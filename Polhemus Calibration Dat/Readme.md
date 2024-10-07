
# Polhemus Calibration and Data Processing Project

## Overview

This project is focused on **calibrating and processing data from Polhemus 3D tracking systems**. It provides a comprehensive set of MATLAB scripts and data files that can:

- Capture motion data in real time
- Process calibration for precise tracking
- Provide useful insights for various experimental setups

Polhemus systems are often used in motion tracking applications such as biomechanics, virtual reality, and 3D spatial mapping. The project implements a custom calibration workflow to optimize the accuracy of the measurements and is divided into distinct modules to handle calibration, data sampling, and other related tasks.

## Directory Structure

The project is organized into various folders and MATLAB scripts. Below is an overview of the structure:

### Main Files

- **polhemusCalData.mat**: Contains core calibration data used during the process.
- **polhemusCalDataSAVE.mat**: Stores backup calibration data.
- **polhemusCalTemp.mat**: Temporary calibration data used for ongoing calculations.
- **SMAT.mat**: Contains additional matrices and variables used for processing.

### `zubs` Folder

The `zubs` folder contains several subdirectories and scripts, each handling different aspects of the calibration and data sampling process:

- **cal**: Contains scripts related to calibration processes.
  - `calFinger.m`: Calibration for finger tracking.
  - `calTable.m`: Main table calibration process.
  - `doELcal.m`: Calibration for EyeLink system integration.
  - `doELdc.m`: Data collection calibration for EyeLink.

- **draw**: Contains scripts for visual output and drawing related to target and feedback displays.
  - `computeRT.m`: Computes real-time feedback.
  - `drawTarg.m`: Draws the target points for calibration.
  - `InstructScreen.m`: Manages screen instructions for the user.
  - `updateFB.m`: Updates feedback on the screen.

- **exps**: Contains experimental scripts related to specific test procedures.
  - `PctVG.m`: Visual guidance percentage test.
  - `Sac.m`: Saccadic experiment script.

- **poll**: Contains scripts for polling and retrieving data from Polhemus or EyeLink systems.
  - `pollPH.m`: Polling data from Polhemus hardware.
  - `pollPHstream.m`: Streaming Polhemus data in real time.

- **samp**: Contains scripts that handle the sampling and saving of data.
  - `datSAVE.m`: Saves sampled data.
  - `sampfixPH.m`: Fixes and processes Polhemus sampled data.

- **ztmp**: Temporary folder to store miscellaneous and intermediate files.
  - `DATmonitor.mat`: Temporary data monitor file.
  - `tmpdat.mat`: Contains intermediate results from real-time calibration.

## Key Features

- **Flexible Calibration**: Custom calibration procedures for Polhemus data that ensure accuracy across various applications.
- **Real-Time Data Polling**: Capability to poll data in real time from Polhemus systems or EyeLink trackers.
- **Modular Structure**: The project is organized to support easy extensions and experimentation with different datasets and setups.
- **Sampling and Feedback**: Provides real-time feedback and stores sampled data during experiments for further analysis.

## Getting Started

To get started with the project, follow these steps:

1. **Clone the Repository**: Clone the repository to your local machine.
   ```
   git clone https://github.com/your-repository-link
   ```

2. **MATLAB Setup**: Ensure you have MATLAB installed with support for the required toolboxes.

3. **Running Calibration**: Use the scripts in the `cal` folder to perform the necessary calibrations based on your Polhemus hardware setup.
   - Example: Running `calTable.m` will calibrate the table-based setup for motion tracking.

4. **Polling Data**: Use the scripts in the `poll` folder to retrieve real-time data.
   - Example: `pollPH.m` polls data from the Polhemus system.
