THE MASTER VERSION OF THIS IS IN THE ABOVE DIRECTORY

===========================
Python code for reading and analyzing DANA data
===========================

LD_Session_Iterator.py - this function iterates through all of the recording session files, runs an analysis function (e.g., Spike_by_Time_Matrices.py), and then stores the output from that analysis into a pkl file in a designated output directory. After running through all sessions, this function can also run the 'meta' higher level analysis that integrates results across datasets.

Spike_by_Time_Matrices.py - this function creates neuron (y axis) by time (x axis) matrices for each dataset. The input is the data directory that 

Plot_Placefields.py - plots placefields for neurons in a given recording session.

===========================
INSTALLATION: I used Anaconda. Used Visual Studio Code for coding and debugging.
===========================
For the future tensorflow stuff, in conda, create an environment (import it) from the yml file. The more advanced things require tf2 and ideally the GPU lib.

Windows setup
See the hardware requirements and software requirements listed above. Read the CUDA® install guide for Windows.

Make sure the installed NVIDIA software packages match the versions listed above. In particular, TensorFlow will not load without the cuDNN64_7.dll file. To use a different version, see the Windows build from source guide.

Add the CUDA, CUPTI, and cuDNN installation directories to the %PATH% environmental variable. For example, if the CUDA Toolkit is installed to C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.1 and cuDNN to C:\tools\cuda, update your %PATH% to match:

SET PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.1\bin;%PATH%
SET PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.1\extras\CUPTI\lib64;%PATH%
SET PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.1\include;%PATH%
SET PATH=C:\tools\cuda\bin;%PATH%

NOTE: The last path depends on where you decide to put that cuDNN folder. It's up to you.

===========================
RUNNING AND EXECUTING CODE:
===========================
1) I select the tf2-GPU (or CPU) environment in Anaconda Navigator
2) I start Visual Studio Code FROM WITHIN ANACONDA NAVIGATOR (i forget to do this and it's a hassle.)
3) In visual studio Code WORKSPACE: select Pyhton 3.7.. tf2-GPU environment/workspace from the pulldown
cutting and pasting in command window???? works but really?

4) on the console, type python (not python3). this will load the approparite python command window.

5) Probably easiest to use Github Desktop to manage updates as I had an issue re-loggin in to Git in VSC.

4) Hitting the run triangle in th upper right is a bit probelmatic as it often seems to like to run Jupyter. I don't typically want Jupyter.
