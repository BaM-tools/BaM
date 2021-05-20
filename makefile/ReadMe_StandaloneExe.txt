Command 'ldd' lists dependencies of an executable, and for BaM the result is:

ldd BaM
	linux-vdso.so.1 (0x00007ffed131e000)
	libgfortran.so.5 => /usr/lib/x86_64-linux-gnu/libgfortran.so.5 (0x00007fcb15442000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fcb152f3000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fcb152d8000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fcb150e6000)
	libquadmath.so.0 => /usr/lib/x86_64-linux-gnu/libquadmath.so.0 (0x00007fcb1509c000)

The first one is linux itself I presume, so we're left with 5 libraries ('shared object' so <=> DLL in Windows, https://stackoverflow.com/questions/9688200/difference-between-shared-objects-so-static-libraries-a-and-dlls-so) that BaM is expecting to find in one of the folders searched at run time (https://unix.stackexchange.com/questions/22926/where-do-executables-look-for-shared-objects-at-runtime). I can see 3 options to avoid problems due to missing .so files:

Option 1: standalone executable: the 5 libraries are incorporated into the executable ("static linking"). However, this requires having the stating verions of each library, i.e. libXXX.a instead of libXXX.so. The list of libraries has to be added at the end of the linker command, i.e. gfortran -o BaM module1.o ... moduleN.o main.o lib1.a lib2.a etc. 
Pros: standalone + supposed to run faster than "dynamic linking" (i.e. looking for .so libraries at run-time as described above).
Cons: larger executable (but no big deal here, just a few Mo) + need to have all static libraries libXXX.a. This is unfortunately not the case here: the 'standalone' executable still expects to find libgcc_s.so.1 (I couldn't find libgcc_s.a, only libgcc.a) and libc.so.6 (I did find a file libc.a, but I guess it's not the good one?). So the 'standalone' executable may or may not work depending on the availability of these two .so libraries on the host machine.

Option 2: provide the 5 libXXX.so files and tell the user to copy the missing ones in one of the paths searched at run time, typically /lib or /usr/lib. Other paths can be defined by the user by modifying the environmental variable LD_LIBRARY_PATH or the content of the file /etc/ld.so.conf
Pros: not our problem any more...
Cons: ugly and geeky from a user's perspective + might be cumbersome in the case or a server where you can't do it by yourself.

Option 3: provide the 5 libXXX.so files in a folder lib_so that stands e.g. next to the executable and let the executable know that .so files are in this folder using the rpath (run-time search path, https://en.wikipedia.org/wiki/Rpath). This can be done with the command: -Wl,-rpath,'$ORIGIN/lib_so'
Pros: solves the cons of option 2.
Cons: if the executable is moved, then folder lib_so should be moved along.



