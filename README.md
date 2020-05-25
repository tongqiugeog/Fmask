# Fmask
The software called Fmask (Function of mask) is used for automated **clouds**, **cloud shadows**, **snow**, and **water** masking for Landsats 4-8 and Sentinel 2 images.

If you have any questions, please contact Zhe Zhu (zhe@uconn.edu) and Shi Qiu (shi.qiu@uconn.edu) at Department of Natural Resources and the Environment, University of Connecticut.

**IMPORTANT:**

This Github page **ONLY** includes the Matlab code for Fmask 4.2. **The Matlab package with GLOBAL AUXILIARY DATA (~1G)** is available at this [Google Drive](https://drive.google.com/drive/folders/1bVwvlGDFOsWnVj5b3MqI5yqRDoi8g935?usp=sharing), where **autoFmask** is the main function for processing an image. **autoFmaskBacth** can process all Landsats 4-8 and Sentinel-2 images into a folder. Note that Mapping Toolbox in Matlab is required for using the source Code.

**Fmask 4.2 softwares** (including standalones with GUI and without GUI) on Windows, Mac, and Linux (~1G without MCR and ~3G with MCR) are also ready to use now! It can be downloaded at this [Google Drive](https://drive.google.com/drive/folders/1bVwvlGDFOsWnVj5b3MqI5yqRDoi8g935?usp=sharing).

This 4.1 version has substantial better cloud, cloud shadow, and snow detection results for Sentinel 2 and better results (compared to the 3.3 version that is currently being used by USGS as the Collection 1 QA Band) for Landsats 4-8 . This one software can be used for **Landsats 4-8 Collection 1 (or Collection 2) Level 1 product (Digital Numbers)** and **Sentinel 2 Level-1C product (Top Of Atmosphere reflectance)** at the same time.


**IMPORTANT:** Majority of the current Collection 1 Landsats 4-8 QA Band provided by USGS are derived form **3.3 Version of Fmask algorithm** based on default parameters (cloud probability is 22.5% and buffer pixel size is 3). For example, (1) The Cloud (bit 4) is based on Fmask cloud mask (0 is not cloud and 1 is cloud in Fmask); (2) The Cloud Confidence (bits 5-6) is based on Fmask cloud probability in which >22.5% is high (11), >12.5% is medium (10), and <12.5% is low (01) with 00 kept for future use; (3) Snow/ice Confidence (bits 9-10) and Cloud Shadow Confidence (bits 7-8) has only low confidence (01) and high confidence (11) which correspond to no and yes respectively in snow/ice and cloud shadow mask provided by Fmask.

**IMPORTANT:** When making the accuracy assessment for Fmask, please dilate 3 pixels for cloud shadow, but no dilation for cloud, snow, and water.


# 4.2 Version
1) Update Fmask tool for processing Landsat Collection 2 data; and allow the cloud probability thershold in a larger range such as [-100, 100] in the GUI version. (Shi Qiu 4/8/2020)

----- 4.1 version below ---- (4.1 version can be download at this [Google Drive](https://drive.google.com/open?id=1l84t_lbp5Cp5v8L_Rbzk3WQHap7MdFwc) for Matlab code and this [Google Drive](https://drive.google.com/drive/folders/1oVefP9G-TD2vhoCaaKCxQjvAnUlrwB19?usp=sharing) for standalones)

2) The cloud shadow mask over water would not be provided at default settings since this will be less meaningful to use and very time-consuming to process. At the same time, fixed the bug that the auxiliary data may not be used for some Sentinel-2 images, of which the extent in the metadata is defined in [0 360] rather than [-180 180]. (Shi Qiu 3/17/2020)

----- 4.0 version below ---- (4.0 version can be download at this [Google Drive](https://drive.google.com/open?id=1SXBnEBDJ1Kbv7IQ9qIgqloYHZfdP6O1O))

3) Fixed the bug that the cloud shadows in Sentinel-2 imagery would be projected along a wrong direction when solar azimuth angle > 180 degrees. (Shi Qiu 01/19/2019)

4) Integrated Cloud Displacement Index (CDI) into this Fmask 4.0 for better seperating clouds from bright surfaces especail for Sentinel-2. The CDI was specially designed to separate clouds from bright surfaces based on the view angle parallax of the three near infrared bands (band 7, 8 and 8a) ([Frantz et al., 2018](https://doi.org/10.1016/j.rse.2018.04.046)).  (Shi Qiu and Zhe Zhu 06/03/2018)

5) Revised the method to identify the potential false positive cloud pixels.  (Shi Qiu and Zhe Zhu 05/23/2018)

6) Restricted the height of the clouds located in the scene boundary into the predicted cloud height derived from its neighboring clouds.  (Shi Qiu 04/05/2018)

7) Removed the overlap between the predicted cloud shadow and the potential cloud shadow layer for cloud shadow detection. (Shi Qiu and Zhe Zhu 03/29/2018)

8) Fixed the bug that the reading blue band using GRIDobj may lead to Nan value for Landsat images. (Shi Qiu 03/26/2018)

9) Improved the computational efficiency specially for cloud shadow matching procedure.  (Zhe Zhu and Shi Qiu 03/24/2018)

10) Released Fmask 4.0 beta version. (Shi Qiu, Zhe Zhu, and Binbin He 03/22/2018)


Please cite the following papers:

**paper 1: Qiu S., Zhu Z., and He B., Fmask 4.0: Improved cloud and cloud shadow detection in Landsats
4-8 and Sentinel-2 imagery, Remote Sensing of Environment, (2019), [doi.org/10.1016/j.rse.2019.05.024](https://doi.org/10.1016/j.rse.2019.05.024) (paper for Fmask 4.0).**

**paper 2: Zhu, Z. and Woodcock, C. E., Improvement and Expansion of the Fmask Algorithm: Cloud, Cloud Shadow, and Snow Detection for Landsats 4-7, 8, and Sentinel 2 images, Remote Sensing of Environment (2014) [doi:10.1016/j.rse.2014.12.014](https://doi:10.1016/j.rse.2014.12.014) (paper for Fmask version 3.2).**

**paper 3: Zhu, Z. and Woodcock, C. E., Object-based cloud and cloud shadow detection in Landsat imagery, Remote Sensing of Environment (2012), [doi:10.1016/j.rse.2011.10.028](https://doi:10.1016/j.rse.2011.10.028) (paper for Fmask version 1.6).**

**paper 4: Qiu S., He B., Zhu Z., et al. Improving Fmask cloud and cloud shadow detection in mountainous area for Landsats 4–8 images, Remote Sensing of Environment (2017), [doi.org/10.1016/j.rse.2017.07.002](https://doi.org/10.1016/j.rse.2017.07.002) (paper for Mountainous Fmask ([MFmask](https://github.com/qsly09/MFmask)), that has been integrated into this Fmask 4.0).**

**paper 5: Qiu, S., Lin Y., Shang R., Zhang J., Ma L., and Zhu Z., Making Landsat Time Series Consistent: Evaluating and Improving Landsat Analysis Ready Data, Remote Sensing (2019), [doi.org/10.3390/rs11010051](https://doi.org/10.3390/rs11010051) (First paper introducing Fmask 4.0 for improving LTS consistency).**


**The training data** for the Fmask 4.0 are availble at the following link: https://landsat.usgs.gov/landsat-8-cloud-cover-assessment-validation-data

**The validation data** for the Fmask 4.0 are comming soon.


After running Fmask 4.0, there will be an image called XXXFmask.tif. The image values are presenting the following classes:

0 => clear land pixel

1 => clear water pixel

2 => cloud shadow

3 => snow

4 => cloud

255 => no observation


# 3.3 Version
Updates (since 3.2):

1) Bug fixed in cloud and cloud shadow matching algorithm.

2) Bug fixed in building 3D cloud objects for small clouds (radium <= 3) for versions where TIRS band is used.

The 3.3 version of **Matlab code** for **Landsats 4-8 in which Landsat 8 has valid TIRS band** can be downloaded at this [link] (https://www.dropbox.com/sh/riruwk721zbl0he/AAAe_ccQiNS7_wHNC3HadOqRa?dl=0)

The 3.3 version of **Windows stand alone software** for **Landsats 4-8 in which Landsat 8 has valid TIRS band** can be downloaded at this [link] (https://www.dropbox.com/sh/ylzub1uzosqidwy/AAC3zmk4M3DSbSoS2OLhR5r9a?dl=0) (provided by Sean Griffin segriffin@gmail.com)

The 3.3 version of **Matlab code** for **Landsats 4-8 in which Landsat 8 does not have valid TIRS band** (zeor values) can be downloaded at this [link] (https://www.dropbox.com/sh/nqepwqp0oo53iio/AAD5abJlu6dt9yg2SIDdfU4Qa?dl=0)

The 3.3 version of **Windows stand alone software** for **Landsats 4-8 in which Landsat 8 does not have valid TIRS band** (zeor values) can be downloaded at this [link] (https://www.dropbox.com/sh/760d18apsfrph3i/AADwp4x0o8XMys6CDBpX6Fgma?dl=0) (provided by Sean Griffin segriffin@gmail.com)

The 3.3 version of **Matlab code** for **Sentinel 2** images can be downloaded at this [link] (https://www.dropbox.com/sh/df5svcpmfddefxu/AACOpdjydXQWmtUDYP1_mDaxa?dl=0) (Fmask 3.3 version desgined for Sentinel-2 data with help from Martin Claverie, Shi Qiu, and Xiaojing Tang). Note that this code is exactly the same as Zhu et al., 2014 proposed. 

Please cite the following papers:

**paper 1: Zhu, Z. and Woodcock, C. E., Object-based cloud and cloud shadow detection in Landsat imagery, Remote Sensing of Environment (2012), doi:10.1016/j.rse.2011.10.028 (paper for Fmask version 1.6.).**

**paper 2: Zhu, Z. and Woodcock, C. E., Improvement and Expansion of the Fmask Algorithm: Cloud, Cloud Shadow, and Snow Detection for Landsats 4-7, 8, and Sentinel 2 images, Remote Sensing of Environment (2014) doi:10.1016/j.rse.2014.12.014 (paper for Fmask version 3.2.).**

The cloud and cloud shadow manual masks used for validating the Fmask mask are availble at the following link:
http://landsat.usgs.gov/ccavds.php

After running Fmask, there will be an image called XXXFmask that can be opened by ENVI. The image values are presenting the following classes:

0 => clear land pixel

1 => clear water pixel

2 => cloud shadow

3 => snow

4 => cloud

255 => no observation

**HOW TO USE**

**Matlab**

Need to install Matlab and have image process and statistics toolboxes and runs on Linux 64 bits machine with 4G+ memory. It can be download and used by the following steps:

1. Download the Matlab code for Fmask 3.3 version by this link and unzip the Fmask folder.

2. Use "addpath" in Matlab environment for the Fmask folder.

3. Type "autoFmask" in the command window.

**Windows Executable**

The instructions can be found in this [link] (https://www.dropbox.com/s/gx7x7htynxk5ulp/Fmask_Windows_Standalone_Instructions_Website.pdf?dl=0). 

**CFmask**
There is also a C version of Fmask performed by USGS. See their site for details [here] (https://github.com/USGS-EROS/espa-cloud-masking).


# 3.2 Version
Updates (since 2.2):

1) Detecting clouds for Landsat 8 using new bands (Zhe 09/04/2013)

2) Remove high probability clouds to reduce commission error (Zhe 09/11/2013)

3) Fix bugs in probs < 0 (Brightness_prob & wTemp_prob) (Zhe 09/11/2013)

4) Add customized snow dilation pixel number (Zhe 09/12/2013)

5) Fix problem in snow detection because of temperature screen (Zhe 09/12/2013)

6) Remove default 3 pixels snow dilation (Zhe 09/20/2013)

7) Fix bug in calculating r_obj and change num_pix value (Zhe 09/27/2013)

8) Remove majority filter (Zhe 10/27/2013)

9) Add dynamic water threshold (Zhe 10/27/2013)

10) Exclude small cloud object < 3 pixels (Zhe 10/27/2013)

**Matlab**
Need to install Matlab and have image process and statistics toolboxes and runs on Linux 64 bits machine with 4G+ memory. It can be download and used by the following steps:

1. Download the Matlab code for Fmask 3.2 version by this link and unzip the Fmask folder.

2. Use "addpath" in Matlab environment for the Fmask folder.

3. Type "autoFmask" in the command window.

**Linux Executable**
Stand alone Linux executable Fmask software which do not need to install Matlab or R and runs on Linux 64 bits machine with 4G+ memory. It is based on the same Fmask 3.2sav Matlab code and it can be download and used by the following steps:

1. Download Fmask 3.2 version Linux package "Fmask_pkg.zip" Use any Brosweer and go to the following [ftp sites]: (http://ftp-earth.bu.edu/public/zhuzhe/Fmask_Linux_3.2v/)

2. Unzip the software using "unzip Fmask_pkg.zip"

3. There will be a new file called MCRInstaller.zip at the same folder and unzip this file.

4. Install MCRInstaller by typing "./install" in the same folder

5. There will be wizard that help you install and there will be two environment variables called "LD_LIBRARY_PATH" and "XAPPLRESDIR" showed up in the wizard. Copy the two variables.

For example, This is what I got:

"On the target computer, append the following to your LD_LIBRARY_PATH environment variable:

/home/amd64

Next, set the XAPPLRESDIR environment variable to the following value:

/home/app-defaults"

6. Edit your .cshrc (.tcsh is the same, for .bash replace it with export LD_LIBRARY_PATH="...") file and add this

"setenv LD_LIBRARY_PATH /home/amd64"

"setenv XAPPLRESDIR /home/app-defaults"

7. Save the shell or bash script and source it;

8. Copy the "Fmask" software to any location you want (for example "/Tools/Fmask");

9. cd into the folder where Landsat bands and .MTL files downloaded and run Fmask by entering "/Tools/Fmask" in the terminals.

There are four important tuning variables that you can play with:

1) "cldpix" is dilated number of pixels for cloud with default values of 3.

2) "sdpix" is dilated number of pixels for cloud shadow with default values of 3.

3) "snpix" is dilated number of pixels for snow with default values of 0.

4) "cldprob" is the cloud probability threshold with default values of 22.5 (range from 0~100). If you want to use default values "/Tools/Fmask" is enough, if you want to customize your own parameters, you can use "/Tools/Fmask cldpix sdpix snpix cldprob", for example "/Tools/Fmask 3 3 0 22.5" in the terminals

**Windows Executable**
Stand alone Linux executable Fmask software which do not need to install Matlab or R and runs on Linux 64 bits machine with 4G+ memory. It is based on the same Fmask 3.2sav Matlab code and it can be download and used by the following steps:

1. Download Fmask 3.2 version Windows package "Fmask_pkg.exe" Use any Brosweer and go to the following [ftp sites]: (http://ftp-earth.bu.edu/public/zhuzhe/Fmask_Windows_3.2v/)

2. Double click "Fmask_pkg.exe" and install it with wizard.

3. There will be a new file called "Fmask.exe" at the same folder and this is your Fmask software

4. Copy the "Fmask.exe" software to any location you want (for example "c:\Tools");

5. cd into the folder where Landsat bands and .MTL files downloaded and run Fmask by entering "c:\Tools\Fmask" in the Command Prompt you can find in the Accessories.

There are four important tuning variables that you can play with:

1) "cldpix" is dilated number of pixels for cloud with default values of 3.

2) "sdpix" is dilated number of pixels for cloud shadow with default values of 3.

3) "snpix" is dilated number of pixels for snow with default values of 0.

4) "cldprob" is the" cloud probability threshold with default values of 22.5 (range from 0~100). If you want to use default values "c:\Tools\Fmask" is enough, if you want to customize your own parameters, you can use “c:\Tools\Fmask cldpix sdpix snpix cldprob", for example “c:\Tools\Fmask 3 3 0 22.5"in the terminals

**CFmask**
There is also a C version of Fmask 3.2 version performed by USGS. See their site for details [here] (https://github.com/USGS-EROS/espa-cloud-masking).

# 2.2 Version
Updates (since 2.1):

1) Fixed a bug find in writing the ENVI header for the UTM zone number (Zhe 02/26/2013).

2) Better cloud and cloud shadow matching results (Zhe 03/01/2013)

Linux Executable
Stand alone Linux executable Fmask software which do not need to install Matlab or R and runs on Linux 64 bits machine with 4G+ memory. It is based on the same Fmask 2.2sav Matlab code and can be downloaded at the "Downloads" tab with help files in the "Wiki" tab.

Windows Executable
Stand alone Windows executable Fmask software created by Sean Griffin (segriffin@gmail.com) which do not need to install Matlab or R and runs on Windows 64 bits machine with 4G+ memory. It is based on the same Fmask 2.2sav Matlab code and can be downloaded at the "Downloads" tab with help files in the "Wiki" tab.

# 2.1 Version
Updates (since 2.0):

1) Process both the new and old "MTL.txt" metadata (Zhe Zhu 10/18/2012)

Matlab
On August 29, 2012, filenames and metadata files associated with the Landsat Level 1 Products has been changed. Modifications are being made to make filenames, metadata fields, and files consistent for all sensors, including upcoming LDCM (Landsat 8) data products. The previous 1.6.3 and 1.6.3sav Matlab version will not be able to work for Landsat data downloaded later than August 29, 2012.

You need to have Matlab software with statistic and image process toolboxes installed on 64 bits machine with 4+G memory. The Fmask 2.1sav software are located at the "Downloads" tab and the help files are in the "Wiki" tab.

Windows Executable
This website also hosts a stand alone Windows executable Fmask software created by Sean Griffin (segriffin@gmail.com) which do not need to install Matlab or R and runs on Windows 64 bits machine with 4G+ memory. It is based on the same Fmask 2.1sav Matlab code and can be at the "Downloads" tab with help files in the "Wiki" tab.

# 2.0 Version
Updates (since 1.6.3):

1) Process TM and ETM+ images with the new "MTL.txt" metadata (Zhe Zhu 09/28/2012)

2) Change the Fmask band name to "Fmask" (Zhe 09/27/2012)

3) Dilate snow by default 3 pixels in 8 connect directions (by Zhe 05/24/2012)

4) Exclude small cloud object <= 9 pixels (by Zhe 03/07/2012)

Matlab
The Fmask 2.0sav software can only process the new Landsat metadata format downloaded after August 29, 2012. It includes the most recent updates since 1.6.3sav as follows:

# 1.6.3 Version
Matlab
The Fmask 1.6.3 & Fmask 1.6.3sav software can only process the old Landsat metadata format downloaded before August 29, 2012.

Two Matlab versions of Fmask: 1.6.3 Version need to use LEDAPS to prepare inputs for Fmask (works with WO.txt and MTL.txt header); 1.6.3 Stand Alone Version (sav) can be run in Matlab environment directly (works with MTL.txt header at the moment). The Fmask 1.6.3v and 1.6.3sav software are located at the "Downloads" tab and the help files are in the "Wiki" tab.

The Fmask1.6.3 version has been validated with the 142 Landsat reference scenes distributed in different continents and the overall accuracy has improved more than 0.06% (overall accuracy of 96.47%) since the 1.6.0 version explained in following paper.

R
There is also a R script translated from Matlab 1.6.3sav version (They will code the most recent Fmask algorithm soon) performed by Joseph Henry (joseph.henry@sydney.edu.au) and Willem Vervoort (willem.vervoort@sydney.edu.au) at Department of Environmental Sciences, Faculty of Agriculture and Environment, The University of Sydney, Australia. See their site for details here.
