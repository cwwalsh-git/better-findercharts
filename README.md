# fast-finderchart
<pre>
Python command-line script to easily create astronomical findercharts with overlaid, customizable data. 
Findercharts and .log file will save into same directory as the finderchart.py file.

Usage (show this message using python finderchart.py -h):
arguments:
  -h, --help            show this help message and exit
  -p [ [ ...]], --pos [ [ ...]]
                        SIMBAD-resolvable object name/coordinate position, or comma-separated list of positions.
  -s [SURVEY [SURVEY ...]], --survey [SURVEY [SURVEY ...]]
                        Survey to retrieve images from, or list of surveys. -h to view all possible surveys.
  -w , --width          Width and height of image in arcminutes
  -m, --no_markers      Boolean to show or hide overlaid markers on image.
  -f , --file_entry     Default pos argument instead takes .csv file path as string. Each row in .csv should contain a
                        SIMBAD-resolvable name or position. Set=2 if .csv has separate ra and dec columns, otherwise
                        set=1.
  -l, --list_surveys    List all surveys instead of running program.
  -q, --quiet           Returns quiet output
  -a, --use_all_surveys
                        Returns images from ALL surveys in survey list.
  -r, --survey_redundant
                         If survey in list fails, will search again with next survey in list.

Main surveys: DSS, SDSSg, Pan-STARRS1/PS1.

argument examples:
  python finderchart.py -p m83, 182 15, m99 -s PS1 -w 3
      =returns 3-arcmin findercharts at M83, ICRS coordinates "182, 15", and M99; uses survey Pan-STARRS1. 
       Note the lack of comma within the coordinate pair (182 15). This is intentional, and be sure to omit this comma to ensure positions are properly parsed.
  python finderchart.py -p example.csv -s DSS -w 2 -f 2
      =returns findercharts with positions from a .csv file in the same working directory. Note that the -f argument changes the -p argument to instead take a filename.
       In this case, -f 2 means that the .csv has coordinate pairs in the first two columns: ie "ra, dec, 123, 12, 50, 13,..." 
       -f 1 would imply that the .csv has positions in the first column only: ie "positions, M83, 182 15, M99,..." 
</pre>
