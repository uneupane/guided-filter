README file for 'Guided Image Filter'.

CONTENTS

main.c          Main processing file
filter.c	Filter
inout.c		Program to read and write PGM file

INPUT FORM
<guidance_image_filename>:   PGM file
<input_filtering_image_filename>:   PGM file
<local window radius>:   r  
<regularization parameter>:  eps
<output_filename>:  PGM file
	
 Note: the size of each local patch is N=(2r+1)^2 except for boundary pixels