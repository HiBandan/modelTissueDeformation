# vertexModel

![alt text](https://github.com/HiBandan/vertexModel-Static/blob/main/logo/vertexModel-3.0.png)


This is a program for simulating emergence of tissue flow and deformation using 2D-vertex model, where experimentally aquired cell shape from the sagittal section of an animal is approximated by polygons. 

Application: in progress

## Note

Testing and contributing is very welcome, especially if you can contribute with new algorithms and features.

## WINDOWS

  ### install gcc/g++ 
  
  download Mingw: https://nuwen.net/mingw.html
  
      run mingw-xx.exe -> this will create a folder MinGw
      
      copy MinGw to C:\
       
      add environment variable: C:\MinGW or C:\MinGW\bin

  ### install eigen:
  
  download eigen: https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
      extract "eigen-3.4.0" folder and rename it to "Eigen"
  
      copy "Eigen" folder to C:\Program Files 
  
  ### install boost:
  
  download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

      extract "boost_1_77_0" folder and rename it to "Boost"
  
      copy "Boost" to C:\Program Files 
  
  ### get the source code
  
  Option 1: Download the package from the Github: https://github.com/HiBandan/modelTissueDeformation/archive/refs/heads/main.zip


  Option 2: Clone from terminal: 
  
    git clone https://github.com/HiBandan/modelTissueDeformation.git

  ### install Code::Blocks: 
  https://code-blocks.en.uptodown.com/windows
    
  ### run 
  
    open (with Code::Blocks): modelTD_WINDOWS.cbp
  
    Projects -> modelTD_WINDOWS.cbp -> add files
        -> modules/cppFiles (all files)
        -> modules/headFiles (all files)
        -> modules/paramFiles (all files)
        
    settings -> compiler -> search directories (compiler) -> add 
  
        -> C:\Program Files\Eigen 
        -> C:\Program Files\Boost 

    project -> 
    
      set program’s arguments: ./modules/paramFiles/const_PARAMETERS.txt
      set program’s arguments: ./modules/paramFiles/vary_PARAMETERS.txt
      build options -> search directories -> compiler -> add -> modules/cppFiles 
      build options -> search directories -> compiler -> add -> modules/headFiles
      build options -> linker settings -> add -> boost_system
      build options -> linker settings -> add -> boost_filesystem
      build options -> linker settings -> add -> libboost_filesystem
      
    run 

## LINUX

  ### install gcc/g++ 
  
       sudo apt update
       sudo apt install build-essential

  ### install eigen:
  
  download eigen: https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
      extract "eigen-3.4.0" folder and rename it to "Eigen"
  
      copy "Eigen" folder to /usr/local/include
  
  ### install boost:
  
  download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

      extract "boost_1_77_0" folder and rename it to "Boost"
  
      copy "Boost" to /usr/local/include

  ### get the source code
  
  Option 1: Download the package from the Github: https://github.com/HiBandan/modelTissueDeformation/archive/refs/heads/main.zip

  Option 2: Clone from terminal: 
  
    git clone https://github.com/HiBandan/modelTissueDeformation.git
    
  ### install Code::Blocks: 
  
  sudo add-apt-repository ppa:damien-moore/codeblocks-stable

  sudo apt update

  sudo apt install codeblocks codeblocks-contrib
 
  ### run 

    open (with code-block): modelTD_LINUX.cbp
    
    Projects -> modelTD_LINUX.cbp -> add files
      -> modules/cppFiles (all files)
      -> modules/headFiles (all files)
      -> modules/paramFiles (all files)
  
    settings -> compiler -> search directories (compiler) -> add 
  
        -> /usr/local/include/Eigen
        -> /usr/local/include/Boost 
       
    project -> 
    
      set program’s arguments: ./modules/paramFiles/const_PARAMETERS.txt
      set program’s arguments: ./modules/paramFiles/vary_PARAMETERS.txt
      build options -> search directories -> compiler -> add -> modules/cppFiles 
      build options -> search directories -> compiler -> add -> modules/headFiles
      build options -> linker settings -> add -> boost_system
      build options -> linker settings -> add -> boost_filesystem
      build options -> linker settings -> add -> libboost_filesystem
  
    run 
