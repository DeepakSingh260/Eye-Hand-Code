# Polhemus Calibration and Data Processing Project

## Overview

This project focuses on calibrating and processing data from Polhemus 3D tracking systems. It provides a comprehensive set of MATLAB scripts for real-time motion data capture, precise calibration, and experimental setups in fields such as biomechanics, virtual reality, and 3D spatial mapping.

## Directory Structure

.
├── polhemusCalData.mat
├── polhemusCalDataSAVE.mat
├── polhemusCalTemp.mat
├── SMAT.mat
└── zubs
├── cal
│ ├── calFinger.m
│ ├── calTable.m
│ ├── doELcal.m
│ └── doELdc.m
├── draw
│ ├── computeRT.m
│ ├── drawTarg.m
│ ├── InstructScreen.m
│ └── updateFB.m
├── exps
│ ├── PctVG.m
│ └── Sac.m
├── poll
│ ├── pollPH.m
│ └── pollPHstream.m
├── samp
│ ├── datSAVE.m
│ └── sampfixPH.m
└── ztmp
├── DATmonitor.mat
└── tmpdat.mat


## Key Components

### Main Data Files

- **polhemusCalData.mat**: Core calibration data
- **polhemusCalDataSAVE.mat**: Backup calibration data
- **polhemusCalTemp.mat**: Temporary calibration data
- **SMAT.mat**: Additional matrices and variables for processing

### Calibration Scripts (zubs/cal)

- **calFinger.m**: Calibrates finger tracking
- **calTable.m**: Main table calibration process
- **doELcal.m**: EyeLink system integration calibration
- **doELdc.m**: Data collection calibration for EyeLink

### Visualization Scripts (zubs/draw)

- **computeRT.m**: Real-time feedback computation
- **drawTarg.m**: Target point drawing for calibration
- **InstructScreen.m**: User instruction screen management
- **updateFB.m**: On-screen feedback updates

### Experimental Scripts (zubs/exps)

- **PctVG.m**: Visual guidance percentage test
- **Sac.m**: Saccadic experiment script

### Data Polling Scripts (zubs/poll)

- **pollPH.m**: Polhemus hardware data polling
- **pollPHstream.m**: Real-time Polhemus data streaming

### Sampling and Saving Scripts (zubs/samp)

- **datSAVE.m**: Sampled data saving
- **sampfixPH.m**: Polhemus sampled data processing and fixing

### Temporary Data (zubs/ztmp)

- **DATmonitor.mat**: Temporary data monitor file
- **tmpdat.mat**: Intermediate real-time calibration results

## Features

- **Custom Calibration**: Tailored procedures for accurate Polhemus data across various applications
- **Real-Time Data Polling**: Capability to poll data in real-time from Polhemus systems and EyeLink trackers
- **Modular Structure**: Organized for easy extensions and experimentation
- **Real-Time Feedback**: Provides immediate feedback during experiments
- **Data Storage**: Stores sampled data for further analysis

## Getting Started

1. Clone the repository:

2. Ensure MATLAB is installed with necessary toolboxes.

3. Run calibration scripts from the `zubs/cal` folder based on your setup.

4. Use polling scripts in `zubs/poll` for real-time data retrieval.

## Usage Example

To calibrate finger tracking:

1. Open MATLAB and navigate to the project directory.
2. Run the following command:
```matlab
Scal = calFinger;
