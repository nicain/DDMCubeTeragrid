Files=\
	MersenneTwister.h \
	DDMCube_Settings.py \
	DDMCube_Master.py \
	DDMCube.pyx \
	DDMCube_Slave.py \
	setup.py \
	/Users/nicain/Library/Python/src/analysisTools/analysisTools.py \
	/Users/nicain/Library/Python/src/pbsTools/pbsTools.py \
	progressMeter.py \
	product.py

Config=-vz

TarName=DDMCubeTeragrid.tar.gz

tar:
	tar $(Config) -cf $(TarName) $(Files)

pushHyak: tar
	scp $(TarName) nicain@hyak.washington.edu:/usr/lusers/nicain/currentProjects/DDMCube

pullHyak:
	rsync -vprzt --stats --progress nicain@hyak.washington.edu:/gscratch/riekesheabrown/nicain/currentProjects/DDMCube/savedResults-Hyak ./
	consolidateResults