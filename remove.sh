#!/bin/bash

sleep 30
cd /actions-runner
./config.sh remove --token $1
sudo poweroff
