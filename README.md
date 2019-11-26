# How to run mcce

## Get the code

### To download the latest code:

git clone https://github.com/GunnerLab/Stable-MCCE.git

### Compile the code:

Compile:

```
make clean
make
```


Add the executable to your path:

```
export PATH={/path/to/mcce/}bin:$PATH
```


**Troubleshooting:**

On Mac OS X, an explicit memory free for tree is not supported. You may have to change this subroutine in file lib/db.c
```C
/* release database memory */
void free_param() {
   tdestroy(param_root, free);
   return;
}
```
to
```C
/* release database memory */
void free_param() {
   return;
}
```

If you see glib.h error:

Install glib:
```sudo apt-get install libglib2.0-dev```

See the command options:
```pkg-config --cflags --libs glib-2.0```

Add the above output to makefile as compiler option (the next is an example):
```-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  -lglib-2.0```


If you see glib.h error

Install glib:
```sudo apt-get install libglib2.0-dev```

See the command options:
```pkg-config --cflags --libs glib-2.0```

Add the above output to makefile as compiler option (the next is an example):
```-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  -lglib-2.0```

## Prepare a working directory and pdb file 

```
mkdir test_lysozyme