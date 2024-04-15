# MotherMachine_NatMicro2024
Scripts and code to analyse bacterial growth in the mother machine used for the Nature Microbiology publication Osbelt et al. 2024. Details are in manuscript.

Information about the included files:
-------------------------------------------------------------------------------------------------

1.
 The Excel-table **"lysis events R1 and R2.xlsx"** lists the mother lineages according to the phenotypes 'filamentous', 'lysis', 'growth arrest' or 'regrowth' <br>
Sheet 1 contains data from both replicates R1 and R2. <br>
Sheet 2 contains data for R1 <br>
Sheet 2 contains data for R2 <br>

2.
In the **CycleData_20231017_4h_R1.mat** and **CycleData_20231025_4h_R2.mat** -files the data of all detected cell division cycles for replicates R1 and R2, respectively, is saved in a table-like manner (struct). A cylce_list-struct contains the following data for every cycle:

**cycle_id** - unique ID of the cell division cycle <br>
**MultipointID** - ID of multipoint (imaged position on the chip), in which the cycle appeared <br>
**trapID** - ID of trap in the multipoint (numerated from left to right), in which the cycle appeared <br>
**times** - list frames, in which the cycle appeared <br>
**begin** - frame, in which the cell of the cycle was born <br>
**end** - frame, in which the cell of the cycle divided, died, or escaped from the trap <br>
**duration** - total number of frames, in which the cycle appeared from 'begin' to 'end' <br>
**length** - list of lengths [µm] of the cell from 'begin' to 'end' <br>
**areas** - list of areas [µm^2] of the cell from 'begin' to 'end' <br>
**XCentroid** - list of X-centroids of the cell from 'begin' to 'end' (location of the cell within the Kymograph in horizontal direction) <br>
**YCentroid** - list of Y-centroids of the cell from 'begin' to 'end' (location of the cell within the Kymograph in vertical direction) <br>
**motherLineageFlag** - criterium if the cell cycle belongs to a mother lineage (cell located permanently at the bottom of the trap) <br>
**parentID** - ID of the parent cell of the cycle. If the cycle starts at frame 1, the 'parentID' entry is empty <br>
**daughtersIDs** - IDs of daughter ID(s) of the cycle. If the cycle ends without a daughter, the 'daughterIDs' entry is empty <br>
**loglength** - list of logarithmized 'length' entries <br>
**time_h** - 'times' entries in hours <br>
**growthrate** - list of growth rates [1/h], calculated by linear regression from the 'loglengths' entries <br>
**timeCenter** - average time, in which the cell division cycle was observed (mean of 'time_h' entry) <br>
**Location** - average location of the cell division cycle (mean of the 'YCentroid' entries) <br>

3.
The **Analyis** folder contains the following files and MATLAB-Scripts: <br>
**Cycle_Generation_Script.mat** - script for the detection of cell division cycles and generation of cycle data, which is stored the CycleData .mat-files <br>
**Plots_Script.mat** - the script is used to read CycleData .mat-files, the  **"lysis events R1 and R2.xlsx"** table, and to generate plots based on these data <br>

The **Colormaps.mat** - file contains a list of colormaps in .mat-format, which was downloaded from Fabio Crameri's resource: <br>
Crameri, F. (2018), Scientific colour maps, Zenodo, doi:10.5281/
zenodo.1243862


  


Instruction for script usage
-------------------------------------------------------------------------------------------------
The **Cycle_Generation_Script.mat** script is organized in the following sections:

1. Selection of working folder, loading of colormaps

2. (OPTIONAL) Deletion of faulty cycles

If faulty cycles (e.g., merging cells or a cell breaking in parts during a cycle) were detected in the course of processing in section 3, the user can interrupt the execution of section 3 and enter the ID of the multipoint which is currently being process. This will lead to the deletion of all cycles for the corresponding multipoint. After correction of the wrongly segmented part of the binary mask (e.g., by using in-built functions in FIJI to add white pixels (+1) to the binary mask or deleting pixels (x 0)), the automated processing of cycle detection in section 3 has to be started again with the same multipoint ID.

3. Automated tracking
   
After the user selects the multipoint, the script loads every kymograph and the corresponding binary mask, one cycle at a time. The binary masks were previously generated in Ilastik using the pixel classification function. The algorithm extracts the data for each segmented object (cell), and links these data according to cycles. After processing the last kymograph in the multipoint, the kymograph together with color-coded overlays is saved as a .jpg image in the 'Analysis' folder and the user is asked to continue with the next multipoint. Each cycle is plotted in the kymograph with a random color, which helps to identify faulty cycles. When all multipoints are processed like that, the data is saved in the **CycleData.mat** - file

4. Growth rate determination

For each cycle longer which is longer than 3 frames, a linear fit is applied on the logarithmized length data to calculate the growth rates


The **Plots_Script.mat** script is organized in the following sections:

1. Load Cycle Data

The cycle data from the files **CycleData_20231017_4h_R1.mat** and **CycleData_20231025_4h_R2.mat** is loaded and pooled into one combined cycle list.

2. Sort data according to multipoints

The data from continuously enumerated multipoints, representing two relicates, is sorted in two categories, whether it belongs to the condition 'WT' or 'TN'

3. Mean growth rate plot

4. Bar plot of mother lineage phenotypes
