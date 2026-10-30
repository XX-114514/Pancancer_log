#!/usr/bin/env bash
# Add an rclone remote for the Cloudflare R2 migration bucket. Run interactively on each server.
# The secret is read without echo and written with printf (a shell builtin), so it never
# appears on a command line, in shell history or in this repository.
set -euo pipefail
RCLONE="${MIG_RCLONE:-rclone}"
NAME="${1:-r2mig}"
unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY
CONF="$("$RCLONE" config file | tail -n 1)"
mkdir -p "$(dirname "$CONF")"; touch "$CONF"; chmod 600 "$CONF"
if grep -q "^\[${NAME}\]" "$CONF"; then
  echo "remote [$NAME] already exists in $CONF; remove it first or pass another name" >&2; exit 1
fi
read -rp "Cloudflare account ID: " ACCOUNT
read -rp "R2 access key ID: " KEY
read -rsp "R2 secret access key: " SECRET; echo
read -rp "Bucket name: " BUCKET
printf '\n[%s]\ntype = s3\nprovider = Cloudflare\naccess_key_id = %s\nsecret_access_key = %s\nendpoint = https://%s.r2.cloudflarestorage.com\nacl = private\nno_check_bucket = true\n' \
  "$NAME" "$KEY" "$SECRET" "$ACCOUNT" >> "$CONF"
unset SECRET
probe="_probe/$(hostname -s)-$(date +%s)"
echo probe | "$RCLONE" rcat "$NAME:$BUCKET/$probe"
"$RCLONE" cat "$NAME:$BUCKET/$probe" >/dev/null
"$RCLONE" deletefile "$NAME:$BUCKET/$probe"
echo "OK: write/read/delete on $NAME:$BUCKET"
echo "Set MIG_REMOTE=$NAME and MIG_BUCKET=$BUCKET in ~/.config/pancancer/migration.local.sh"
