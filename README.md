# MotherMachine_NatMicro2024
Scripts and code to analyse bacterial growth in the mother machine used for the Nature Microbiology publication Osbelt et al. 2024. Details are in manuscript.

Information about the included files:
-------------------------------------------------------------------------------------------------


The Excel-table "lysis events R1 and R2.xlsx" lists the mother lineages according to the phenotypes 'filamentous', 'lysis', 'growth arrest' or 'regrowth'
Sheet 1 contains data from both replicates R1 and R2.
Sheet 2 contains data for R1
Sheet 2 contains data for R2

In the CycleData_20231017_4h_R1.mat and CycleData_20231025_4h_R2.mat-files the data of all detected cell division cycles for replicates R1 and R2, respectively, is saved in a table-like manner (struct). A cylce_list-struct contains the following data for every cycle:

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


/-----------------------------------------------------------------------------------------------/

bla
-------------------------------------------------------------------------------------------------
