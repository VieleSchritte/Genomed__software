#!/bin/bash

ROOT=/var/www/cog.doccard.club

NAME="paternity_app"                            #Name of the application (*)
DJANGODIR=${ROOT}/${NAME}             # Django project directory (*)
SOCKFILE=${ROOT}/run/gunicorn.sock        # we will communicate using this unix socket (*)
USER=www-data                                        # the user to run as (*)
GROUP=www-data                                     # the group to run as (*)
NUM_WORKERS=2                                    # how many worker processes should Gunicorn spawn (*)
DJANGO_SETTINGS_MODULE=${NAME}.settings             # which settings file should Django use (*)
DJANGO_WSGI_MODULE=${NAME}.wsgi                     # WSGI module name (*)

echo "Starting $NAME as `whoami`"

# Activate the virtual environment
#cd $DJANGODIR
#source /var/www/test/venv/bin/activate
#export DJANGO_SETTINGS_MODULE=$DJANGO_SETTINGS_MODULE
#export PYTHONPATH=$DJANGODIR:$PYTHO:qNPATH

# Create the run directory if it doesn't exist
RUNDIR=$(dirname $SOCKFILE)
test -d $RUNDIR || mkdir -p $RUNDIR

# Start your Django Unicorn
# Programs meant to be run under supervisor should not daemonize themselves (do not use --daemon)
exec /usr/local/bin/gunicorn ${DJANGO_WSGI_MODULE}:application \
  --name $NAME \
  --workers $NUM_WORKERS \
  --user $USER \
  --bind=unix:$SOCKFILE
