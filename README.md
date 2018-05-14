# RayTimes_through3D

USAGE:
./xumich\_cj\_ttc lines raypath input >> out

lines : integer, number of lines in raypath file

raypath : file, four columns output of taup\_path with geographic coorinates.
 
The four columns are:
Epicentral degrees from source, Depth, ray latitude, ray longitude.

input: See input.TEMPLATE for example


OUTPUT:
Three floats fomatted as such:

(F10.3, F10.3, F10.3)                                                            

Theoretical ttime, Tomographic ttime, ttime difference

./do\_compile should produce executable

When installing on other systems, change line 150 and line 441 of 
/src/umich\_cj\_ttc.f90 to the path of the models directory.
