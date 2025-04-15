## MCCSTA
### The codes employed in the 3-D multimodal atlas project of the marmoset cerebral cortex.
### The data can be interactively view at https://db.cngb.org/stomics/mbcsta/.
*We sincerely thank the support in database building and maintaining provided by China National GeneBank (CNGB).*
### The paper can be accessed at ...

*The cell types and connectivity of primate brains remain unclear. We created a 3-D multimodal atlas of the marmoset brain by combining single-cell spatial transcriptomes, MRI, and retrograde labeling.*
*This multimodal atlas offers new understanding of primate cortical organization and evolution.*

![image](https://github.com/user-attachments/assets/a72a2497-bdc6-43ca-bf5c-9255816f3460)

<img src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSxuMOcdNwWiWNb_FPaD-Z54i0oa_R8Qu2H0A&s" alt="Example Image" width="30" height="30">
Code registrating MRI to spatial transcriptome can be found in <u>Registration of Stereo-seq Data with MRI Data</u> directory.

### The process of Stereo-seq analysis of this paper in Stereo_anlysis:

**01.Atlas_plot:** The code that showed the plotting of Stereo-seq cellbin data relates to Figure 1.

**02.Segment2flatmap:** The code used to convert the segment data into segment-flatmap and to combine the cell-type and gene data within the segment-flatmap relates to Figure S1.

**03.Pr-Al_calculation:**  The code that showed the calculation of Pr - Al Index of different data (four species, region - level, segment - level) relates to Figure 2.

**...**

**10.AU_region_gene:**  The code demonstrates the process of showing that the A1 region in humans is more similar to that in marmosets to Figure 7.

**additional.1.MRIflatmaptoshp:**  The code was used to convert the four species MRI gii format (tif format) flatmap into shape file.

**additional.2.MRIflatmapvisual:**  The code that showed how to plot the species shape file of flatmap.
