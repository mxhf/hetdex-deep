#!/bin/bash
ffmpeg -framerate 10 -i rendering%02d.png -c:v libx264 -r 30 -pix_fmt yuv420p -b 6400k out.mp4
