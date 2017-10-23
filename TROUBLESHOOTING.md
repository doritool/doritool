If DoriTool does not run correctly, and give errors, I recommind you to install all again following the next steps and advices.

1. From DoriTool website: https://doritool.github.io/   (setup).
I recommind you to download it using linux.
If you want to use it with windows,  it will work only with windows 10 versions.
(It is necessary to download docker in both cases)
Once you are downloading docker tool, if ask you to download docker toolbox, means that your windows version is not the one required and will not work doritool correctly (we realized that).

2. Once docker it is installed, you must Activate BIOS VT-X/AMD-v if is not enabled. Check it.

3. Download the DoriTool repository directly

  - git clone https://github.com/doritool/doritool.git (with linux)

    or the zip file (with windows)

  - https://github.com/doritool/doritool/archive/master.zip

Notice:

For window’s users download the DoriTool repository using cmd (click windows R and write cmd, the shell will be opened). Or directly from the docker toolbox Shell.

4. Doritool works better using Ethernet than WIFI. (VEP works better in this sense)

5. The input list must be inside DoriTool folder (the one downloaded before). Path to doritool folder using cmd or terminal in case of linux.

6. Using the terminal write:

    - Linux users

    `./doritool -i input_file.txt`

    - Windows users

    `doritool.bat -i input_file.txt`

add the arguments you want to obtain LD, GTEX in the same command line aswell
example:

- Linux users

  `./doritool -l 0.90 -i imputfile.txt`

- Windows users:

  `doritool.bat -l 0.90 -i imputfile.txt`

