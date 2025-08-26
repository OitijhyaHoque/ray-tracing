#!/usr/bin/bash

g++ ./src/2005049_main.cpp ./src/2005049_obj.cpp -o r -lglut -lGLU -lGL && ./r && rm r
