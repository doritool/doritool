If DoriTool does not run correctly, and give errors, I recommind you to install all again following the next steps and advices.

1. From DoriTool website: https://doritool.github.io/   (setup).
I recommend you to download it using any Linux distribution.
If you want to use it with windows, it will work only with windows 10 versions.
(It is necessary to download docker in both cases)
Once you are downloading docker tool, if ask you to download docker toolbox, means that your windows version is not the one required and will not work doritool correctly (we realized that).

2. Once docker it is installed, you must Activate BIOS VT-X/AMD-v if is not enabled. Check it.

Notice:

For window’s users download the DoriTool repository using cmd (click windows R and write cmd, the shell will be opened). Or directly from the docker toolbox Shell.

Doritool works better using Ethernet than WIFI. (VEP works better in this sense)

3. Download/Update the Docker image from Docker Hub with the next command (<font color="red"> the image will be downloaded, be patient</font>)

    `docker pull doritool/doritool`

4. Using the terminal write:

    - Linux users

    `docker run --rm -u "$(id -u)" -v "$(pwd)":/home/vep/doritool doritool/doritool" -i input_file.txt`

    - Windows users

    `docker run --rm -u 1000 -v %cd%:/home/vep/doritool doritool/doritool -i input_file.txt`

add the arguments you want to obtain LD, GTEX in the same command line aswell
example:

- Linux users

  `docker run --rm -u "$(id -u)" -v "$(pwd)":/home/vep/doritool doritool/doritool" -l 0.90 -i input_file.txt `

- Windows users:

  `docker run --rm -u 1000 -v %cd%:/home/vep/doritool doritool/doritool -l 0.90 -i inputfile.txt`

