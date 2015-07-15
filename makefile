all:
	current_dir = $(notdir $(shell pwd)) 
	cd $current_dir/third-party/oases; make 'VELVET_DIR=$current_dir/third-party/'
	cd $current_dir/third-party/bamUtil_1.0.12; make
	cd $current_dir/third-party/bedtools2; make
	cd $current_dir/third-party/khmer; make #not sure if this is required
	cd $current_dir/third-party/subread-1.4.6/src; make -f Makefile.Lnux
	#Not sure what is required for trim galore yet
