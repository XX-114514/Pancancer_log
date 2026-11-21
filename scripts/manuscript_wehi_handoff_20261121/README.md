# Manuscript receiving verifier

This standard-library verifier takes an external handoff root. It validates the archive against its transfer manifest, rejects traversal and symlink ZIP members, extracts into a new directory, verifies every frozen manifest file, and runs the delivered plotting-input hash verifier. Existing files must pass the same checks; they are not overwritten.

```bash
python3 verify_and_prepare.py --handoff-root "${PRIVATE_HANDOFF_ROOT}"
```

Run the receiver in its intended server environment. A local or synthetic execution does not prove remote receipt. The surrounding private transfer wrapper uses existing user-owned R2 configuration and immutable object copies. The manuscript, plotting arrays, configuration and bucket credentials are supplied outside this repository. No plotting/model fit is executed by the receiving verifier.
