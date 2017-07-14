The models are attached in two versions, MAtlab 2007a and MAtlab 2010a.
The difference is in subsystem Display Results.
The 2010a version contains better function "wrapToPi". This function deals with cases when yaw angle 
exceeds values of Pi. This m-function unfortunately can not be implemented in version 2007a, so little delay
causes the four vertical lines touching the value of two Pi. This is not error.

The models were created in version 2007a. So, when you want to run version 2010a, prior to running you have
to save model as version 2010a.

If you have problems with encoding, run following command:

set_param(0, 'CharacterEncoding', 'Windows 1252');

If you have problems with some matrix operations (Simulink blocks from blockset Matrix Operations are not included 
in every version of Simulink) you have to replace these blocks by m-scripts (LU inverse matrix) or by constants (unit matrix).

