Middlelayer
-----------

Installation
- build specmx
- in output move content of +specmx to +specmx/+internal
- move all files from here to +specmx

Comments
--------

- read only properties will throw a MATLAB:class:noSetMethod error when someone tries to set them

- the properties are not cached, everytime they are read, a new call into specmx is made,
  including some performance overhead

- modifications in the constructors of the classes (like default arguments) can lead to unwanted
  behavior because SWIG uses them two times during construction, the second time with a SwigRef
  instance, therefore we need to check we are not in the second call when processing default arguments
  see File for an example how to do that

- in Stack all Properties with vector output give row vectors (underlying functions give column vectors)
  (this is more in line with standard Matlab behavior of size, ...)

- in Stack all Properties with vector input can process row or column vectors (underlying functions
  only accept column vectors)

- Matlab uses delete(handle class instance) to delete the handle (e.g. close files, ...). This is also
  the preferred way here to close a file before the variable holding it falls out of scope