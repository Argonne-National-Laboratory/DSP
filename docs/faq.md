# Frequently Asked Questions

Q: I have a linking error with `libma27.a`, saying "recompile with -fPIC".

A: Please configure `libma27.a` with `-fPIC` flag and compile it again. For example,
```
./configure FFLAGS=-fPIC
make
make install
```
