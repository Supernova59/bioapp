#!/bin/bash
cd ~/bioapp && source .venv/bin/activate
D="${XDG_RUNTIME_DIR:-/tmp}"
if ! pgrep socat > /dev/null; then
  rm -f "$D/discord-ipc-0"
  socat UNIX-LISTEN:"$D/discord-ipc-0",fork EXEC:"/mnt/c/Users/nonol/npiperelay.exe -ep -s //./pipe/discord-ipc-0" &
  sleep 1
fi
python main.py
