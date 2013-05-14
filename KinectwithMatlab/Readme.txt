Release note V1.0


System Requirement:

Windows 7 (x86 or x64)
Microsoft .NET Framework 4.0
Microsoft Visual C++ 2010 compiler
Kinect for Windows SDK beta
Matlab 2010b or above(Which can recognize Microsoft Visual C++ 2010 compiler)


Setting up the environment.

Step 1: 

Go to Kinect SDK download page. Download Kinect for Windows SDK beta with your x64 or x86 computer:
http://research.microsoft.com/en-us/um/redmond/projects/kinectsdk/download.aspx
(Note: My code is developed under Ver 1.0012)

Step 2:

Set up your mex compiler as "Microsoft Visual C++ 2010" in Matlab.
In matlab, type
mex -setup

Step 3:

Change current directory path to the code's directory.
Compile the getimagedata.cpp, specify the include path and the lib file:
I use the typical path as example:
mex getimagedata.cpp -I'C:\Program Files (x86)\Microsoft Research KinectSDK\inc' 'C:\Program Files (x86)\Microsoft Research KinectSDK\lib\MSRKinectNUI.lib'

Usage:

In Matlab console:
type "[a b] = getimagedata(1);" at first time to initialize Kinect.
Then type
while(1); [a b] = Kinect(); subplot(1,2,1); imagesc(a); axis image; subplot(1,2,2); imagesc(b/255); axis image; drawnow; end;

Enjoy!


Bugs already known:

Can't compile more than once in the matlab. Need to restart Matlab.
The shutdown function does not work normally, hope it can be solved in next release of Kinect for Windows SDK beta.
