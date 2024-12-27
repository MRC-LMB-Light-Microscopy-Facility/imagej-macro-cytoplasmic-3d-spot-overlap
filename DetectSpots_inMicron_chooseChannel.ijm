#@Integer(label="Channel 1", value=2) SPOT1
#@Integer(label="Channel 2", value=3) SPOT2
#@Integer(label="Nuclei channel", value=1) NUC
#@Integer(label="Minimum spot size", value=1) minimumSpotSize
#@Integer(label="Maximum spot size", value=1) maximumSpotSize
#@Float(label="Dilation factor",value=0) DilationFactor

/*
 * This script will find spots in 2 channels, count and calculate their overlap.
 * The number of spots is then normalized by the cytoplasmic volume (ie not nuclear)
 * 
 * Requirments: 
 *   image science 
 *   
 * Script written by Ulrike Schulze for Natalia and Feline in 2022.
 */

run("FeatureJ Options", "isotropic progress");
run("3D OC Options", "volume nb_of_obj._voxels dots_size=1 font_size=10 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
getDimensions(width, height, channels, slices, frames);
name = getTitle();

NoOfSpots		=newArray();
SizeFilterMin	=newArray();
SizeFilterMax	=newArray();
NoOfVoxels		=newArray();
CellVolume		=newArray();
SPOTS = newArray(SPOT1, SPOT2);


setBatchMode("hide");

for (spotc = 0; spotc <= 1; spotc++){
	
	ch = SPOTS[spotc];
	
	//Laplace Filter for each zslice
	selectWindow(name);
	run("Duplicate...", "duplicate channels="+ch);
	rename("spot"+ch);
	run("Enhance Contrast", "saturated=0.35");
	
	for (i = 1; i <= slices; i++){
		selectWindow("spot"+ch);
		setSlice(i);
		run("Duplicate...", " ");
		rename(i);
		run("FeatureJ Laplacian", "compute smoothing=0.1");
		close("i");
		if (i>1){
			run("Concatenate...", "  title=[1 Laplacian] open image1=[1 Laplacian] image2=["+i+" Laplacian] image3=[-- None --]");	
		}
		close(i);
	}
	//convert to mask
	setAutoThreshold("Triangle stack");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Triangle background=Light black");
	rename("Mask "+spotc);
	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	run("3D Objects Counter", "threshold=128 slice=12 min.="+minimumSpotSize+" max.="+maximumSpotSize+"  objects statistics summary");
	
	//make sanity check stack
	selectWindow("Objects map of Mask "+spotc);
	setThreshold(1, 65535, "raw");
	run("Convert to Mask", "method=Triangle background=Dark black");
	run("Duplicate...", "duplicate");
	rename("sanity check");
	run("Find Edges", "stack");
	run("16-bit");
	run("Merge Channels...", "c1=[sanity check] c4=[spot"+ch+"] create ignore");
	rename("sanity check Channel"+ch);
	close("Mask "+spotc);
	
	//extract no of voxels	
	selectWindow("Statistics for Mask "+spotc);
	voxels = Table.getColumn("Nb of obj. voxels");
	Array.getStatistics(voxels, minimum, maximum, mean);
	number_of_elements = lengthOf(voxels);
	NoOfVoxels = Array.concat(NoOfVoxels,mean*number_of_elements);
	close("Statistics for Mask "+spotc);	
	logString = getInfo( "Log" );
	
	//extract no of spots
	a =  indexOf(logString, ": ")+1;
	b =  indexOf(logString, "objects");
	c= substring(logString,a+1,b);
	c = parseInt(c);
	NoOfSpots=Array.concat(NoOfSpots,c);
	close("Log");

}

// Calculate overlap
selectWindow("Objects map of Mask 0");
for (i = 1; i <= DilationFactor; i++){
	run("Erode", "stack");
}
imageCalculator("AND create stack", "Objects map of Mask 0","Objects map of Mask 1");
rename("AND");

getStatistics(area, mean, min, max, std, histogram);
if (max==0) {
	run("3D Objects Counter", "threshold=0 slice=12 min.=1 max.=1E6 exclude_objects_on_edges objects statistics summary");} //if there's no overlap
else {
	run("3D Objects Counter", "threshold=128 slice=12 min.=1 max.=1E6  objects statistics summary");//if there's  overlap
}
	
// extract no of voxels	
selectWindow("Statistics for AND");	
voxels = Table.getColumn("Nb of obj. voxels");
Array.getStatistics(voxels, minimum, maximum, mean);
number_of_elements = lengthOf(voxels);
NoOfVoxels = Array.concat(NoOfVoxels,mean*number_of_elements);
close("Statistics for AND");	

logString = getInfo( "Log" );
// extract no of spots
a =  indexOf(logString, ": ")+1;
b =  indexOf(logString, "objects");
c = substring(logString,a+1,b);
c = parseInt(c);
NoOfSpots=Array.concat(NoOfSpots,c);
close("Log");

close("Objects map of Mask 0");	
close("Objects map of Mask 1");	
close("Objects map of AND");
close("Log");

selectWindow("AND");
rename("overlap");	

//normalise by cell volume (==everything that is not a nucleus)
selectWindow(name);
run("Duplicate...", "duplicate channels="+NUC);
rename("Channel 4");
//if  2d image, the thresholding works slightly different, has a different name in the end
if (slices==1){
	run("Duplicate...", "title=[MASK_Channel 4]");
	setAutoThreshold("Otsu dark stack");
	run("Convert to Mask", "method=Otsu background=Dark create");
} else {
	setAutoThreshold("Otsu dark stack");
	run("Convert to Mask", "method=Otsu background=Dark create");
}
run("Fill Holes", "stack");
run("Invert", "stack");

run("3D Objects Counter", "threshold=128 slice=12 min.=1 max.=1E60 statistics summary");
selectWindow("Statistics for MASK_Channel 4");
Volume = Table.getColumn("Volume (micron^3)");
Array.getStatistics(Volume, minimum, maximum, mean);
number_of_elements = lengthOf(Volume);
CellVolume = Array.concat(CellVolume,mean*number_of_elements);
close("Statistics for MASK_Channel 4");	
selectWindow("MASK_Channel 4");
run("Find Edges", "stack");
run("16-bit");
run("Merge Channels...", "c1=[MASK_Channel 4] c4=[Channel 4] create ignore");
rename("sanity check Channel 4");
close("Log");	




//make Result table
if (isOpen("Results") ==0){
	Table.create("Results")};

a = nResults;


setResult("No of Spots Ch"+SPOTS[0]+"/ micron^3", a, NoOfSpots[0]/CellVolume[0]);
setResult("No of Spots Ch"+SPOTS[1]+"/ micron^3", a, NoOfSpots[1]/CellVolume[0]);
setResult("No of Spots Overlap/ micron^3", a, NoOfSpots[2]/CellVolume[0]);
setResult("Manders M1: Overlap/Ch"+SPOTS[0], a, NoOfVoxels[2]/NoOfVoxels[0]);
setResult("Manders M2: Overlap/Ch"+SPOTS[1], a, NoOfVoxels[2]/NoOfVoxels[1]);
setResult("No of Spots Channel "+SPOTS[0], a, NoOfSpots[0]);
setResult("No of Spots Channel "+SPOTS[1], a, NoOfSpots[1]);
setResult("No of Spots Channel Overlap", a, NoOfSpots[2]);
setResult("Cell Volume in micron^3", a, CellVolume[0]);



setResult("Image", a, name);
setResult("Size Filter Min", a, minimumSpotSize);
setResult("Size Filter Max", a, maximumSpotSize);
setResult("Dilation ", a, DilationFactor);

setBatchMode("exit and display");

run("Tile");
run("Synchronize Windows");

