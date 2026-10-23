"""Read one user-provided file into RAM, with complete identity before parsing.

The caller supplies the frozen path and public whole-file identity. This module
does not download, write, rename, chmod or parse the input. Short or changing
uploads never produce an accepted buffer.
"""
import hashlib
import io
import os
import stat
import time


class UploadedInputError(ValueError):
    def __init__(self, code_name):
        super().__init__(code_name)
        self.code_name = code_name


def signature(info):
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns,
            info.st_ctime_ns, info.st_uid, info.st_mode)


def read_verified_upload(path, *, expected_bytes, expected_MD5, expected_SHA256,
                         maximum_input_bytes, maximum_seconds, progress=None,
                         chunk_bytes=1024**2):
    if not (0 < expected_bytes <= maximum_input_bytes and chunk_bytes > 0
            and maximum_seconds > 0):
        raise UploadedInputError('INVALID_INPUT_RESOURCE_BUDGET')
    start = time.monotonic()
    before = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(before.st_mode) or before.st_uid != os.geteuid():
        raise UploadedInputError('INPUT_NOT_OWNED_REGULAR_FILE')
    if before.st_size != expected_bytes:
        raise UploadedInputError('UPLOAD_LENGTH_NOT_COMPLETE')
    flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
    descriptor = os.open(path, flags)
    buffer = io.BytesIO()
    try:
        opened = os.fstat(descriptor)
        if signature(opened) != signature(before):
            raise UploadedInputError('INPUT_CHANGED_BEFORE_READ')
        md5, sha = hashlib.md5(), hashlib.sha256()
        count = 0
        with os.fdopen(descriptor, 'rb', closefd=False) as stream:
            while count < expected_bytes:
                if time.monotonic() - start > maximum_seconds:
                    raise UploadedInputError('INPUT_READ_TIME_BUDGET')
                chunk = stream.read(min(chunk_bytes, expected_bytes-count))
                if not chunk:
                    raise UploadedInputError('INPUT_TRUNCATED_DURING_READ')
                buffer.write(chunk); md5.update(chunk); sha.update(chunk)
                count += len(chunk)
                if progress is not None:
                    progress(count)
            if stream.read(1):
                raise UploadedInputError('INPUT_GREW_DURING_READ')
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
        if signature(opened) != signature(after) or signature(opened) != signature(current):
            raise UploadedInputError('INPUT_CHANGED_DURING_READ')
        if md5.hexdigest() != expected_MD5 or sha.hexdigest() != expected_SHA256:
            raise UploadedInputError('WHOLE_FILE_IDENTITY_MISMATCH')
        # Buffered reads and coarse/cached filesystem timestamps can hide a
        # rewrite of a prefix already copied to RAM. Verify the same descriptor
        # again without allocating another file-sized buffer.
        os.lseek(descriptor, 0, os.SEEK_SET)
        repeated_md5, repeated_sha = hashlib.md5(), hashlib.sha256()
        repeated_count = 0
        while repeated_count < expected_bytes:
            if time.monotonic() - start > maximum_seconds:
                raise UploadedInputError('INPUT_READ_TIME_BUDGET')
            chunk = os.read(descriptor, min(chunk_bytes, expected_bytes-repeated_count))
            if not chunk:
                raise UploadedInputError('INPUT_CHANGED_DURING_READ')
            repeated_count += len(chunk)
            repeated_md5.update(chunk); repeated_sha.update(chunk)
        if (os.read(descriptor, 1) or repeated_md5.hexdigest() != expected_MD5
                or repeated_sha.hexdigest() != expected_SHA256
                or signature(os.fstat(descriptor)) != signature(opened)
                or signature(os.stat(path, follow_symlinks=False)) != signature(opened)):
            raise UploadedInputError('INPUT_CHANGED_DURING_READ')
        if time.monotonic() - start > maximum_seconds:
            raise UploadedInputError('INPUT_READ_TIME_BUDGET')
        proof = dict(status='COMPLETE_USER_UPLOAD_READ_ONLY_RAM_INPUT_VERIFIED',
            exact_length_before_and_after=True, owned_regular_file=True,
            same_descriptor_and_path_inode=True, unchanged_size_mtime_ctime_and_mode=True,
            whole_official_MD5_verified=True, whole_frozen_SHA256_verified=True,
            second_descriptor_pass_whole_hash_verified=True,
            no_network_requests=True, no_input_write_rename_or_permission_change=True,
            no_additional_native_or_expanded_payload_files=True,
            bytes_read=count, existing_user_file_retained=True)
        return buffer, md5, sha, proof
    except BaseException:
        buffer.close()
        raise
    finally:
        os.close(descriptor)
