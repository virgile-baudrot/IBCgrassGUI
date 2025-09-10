#!/usr/bin/env bash
set -e

# --- Xvfb + Xfce + VNC/noVNC ---
export DISPLAY=:1
export XDG_RUNTIME_DIR=/tmp/runtime-root
mkdir -p "$XDG_RUNTIME_DIR" && chmod 700 "$XDG_RUNTIME_DIR"

Xvfb :1 -screen 0 1280x800x24 +extension RANDR &
sleep 1
dbus-launch startxfce4 >/tmp/xfce.log 2>&1 || true &

x11vnc -display :1 -forever -nopw -shared -rfbport 5900 -noxdamage >/tmp/x11vnc.log 2>&1 &
websockify --web=/usr/share/novnc/ 6080 localhost:5900 >/tmp/websockify.log 2>&1 &

# --- App ---
APP_DIR="/opt/IBCgrassGUI"
echo "[start.sh] Using APP_DIR=$APP_DIR" | tee /tmp/ibcgrass.log
ls -la "$APP_DIR" | tee -a /tmp/ibcgrass.log

# Normalize CRLF -> LF (since copy from Windows)
if command -v dos2unix >/dev/null 2>&1; then
  find "$APP_DIR" -maxdepth 1 -type f \( -name "RunIBCgrassGUI_Linux.sh" -o -name "RunIBCwithoutGUI.R" \) -print0 \
    | xargs -0 dos2unix -k >/dev/null 2>&1 || true
  find "$APP_DIR/Model-files" -type f -name "*.sh" -print0 2>/dev/null \
    | xargs -0 dos2unix -k >/dev/null 2>&1 || true
else
  sed -i 's/\r$//' "$APP_DIR/RunIBCgrassGUI_Linux.sh" 2>/dev/null || true
  sed -i 's/\r$//' "$APP_DIR/RunIBCwithoutGUI.R" 2>/dev/null || true
fi

# IMPORTANT: remove to force new compilation
rm -f "$APP_DIR/Model-files/IBCgrassGUI" 2>/dev/null || true

# Launch script Linux
cd "$APP_DIR"
if [ -f ./RunIBCgrassGUI_Linux.sh ]; then
  chmod +x ./RunIBCgrassGUI_Linux.sh || true
  echo "[start.sh] Launching RunIBCgrassGUI_Linux.sh (will compile under ./Model-files/)" | tee -a /tmp/ibcgrass.log
  bash ./RunIBCgrassGUI_Linux.sh >>/tmp/ibcgrass.log 2>&1 &
else
  echo "WARNING: Script missing: RunIBCgrassGUI_Linux.sh" | tee -a /tmp/ibcgrass.log
fi

# Keep Container alive + print logs
tail -F /tmp/websockify.log /tmp/x11vnc.log /tmp/ibcgrass.log /tmp/xfce.log
