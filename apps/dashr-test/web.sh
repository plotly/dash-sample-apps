#!/bin/sh
if [ "$LANGUAGE" = "PYTHON" ]
then gunicorn app:runserver
else fakechroot fakeroot chroot /app/.root /bin/sh -c 'cd /app && /usr/bin/R -f /app/run.R --gui-none --no-save'
fi
