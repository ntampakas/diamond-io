#!/bin/bash

sleep 30
cd /actions-runner
./config remove --token $1
sudo poweroff
