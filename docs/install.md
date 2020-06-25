# How to Install

This describes the following:

* how to correctly clone the repository with submodules
* the prerequisites for installation
* how to install DSP and Julia interface

## Download

You can clone this repository in your preferred directory by typing

```shell
git clone --recursive https://github.com/Argonne-National-Laboratory/DSP.git
```

or

```shell
git clone https://github.com/Argonne-National-Laboratory/DSP.git
cd DSP
git submodule update --init --recursive
```

## Build and Setting

Please set `UserConfig.cmake` as follows.
Assuming that you are at the root directory of DSP, type

```shell
mkdir build
cd build
cmake ..
```

A shared object is installed in ``./lib`` directory. 
Once the installation has been successfully done, you need to set environment variable ``(DY)LD_LIBRARY_PATH``.
