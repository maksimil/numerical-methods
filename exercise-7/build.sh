#!/usr/bin/env sh

logstatus() {
    echo "\e[1;34m>> $1\e[0m"
}

command=$1
config=$2

if [ "$command" != "build" ] && [ "$command" != "run" ]; then
    logstatus "Unknown command"
    exit 1
fi

if [ "$command" = "build" ] || [ "$command" = "run" ]; then
    logstatus "Building the app"
    if [ "$config" = "release" ]; then
        cmake -Bbuild . -DCMAKE_BUILD_TYPE=Release
    elif [ "$config" = "debug" ]; then
        cmake -Bbuild . -DCMAKE_BUILD_TYPE=Debug
    else
        logstatus "Unknown build config"
        exit 1
    fi
    make -Cbuild
fi


if [ "$command" = "run" ]; then
    logstatus "Running the app"
    ./build/gtktest
    logstatus "App exited with $?"
fi
