#!/bin/bash

cd doc && make html && make latexpdf && cd ..

PYTHONPATH="`pwd`:`pwd`/modules/:$PYTHONPATH"
NOSETESTS=`which nosetests`

if [[ ! -f "$NOSETESTS" ]] ; then
    NOSETESTS=`which nosetests2`
fi

if [[ ! -f "$NOSETESTS" ]] ; then
    echo 'Cannot find nosetests or nosetests2';
else
   echo "Using $NOSETESTS"
   $NOSETESTS --with-doctest --with-coverage \
              --exe \
              --cover-package="etrsitrs" \
              --cover-tests \
              --cover-html \
              --cover-html-dir=coverage \
              --cover-erase \
              -x $@
fi
