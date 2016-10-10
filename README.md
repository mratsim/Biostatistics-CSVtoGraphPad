# Biostatistics-CSVtoGraphPad
From hundreds of flow cytometry CSV file to cleaned up Excel file, ready for GraphPad copy-paste.


This was used to help my wife manage all the data she had during her PhD. (more than 130 days of experiments, several GB of data files).

Everything was stored in CSV file, challenge was extracting relevant information (date, cell line, ...) from file name and CSV row.
Plus there were a couple error in label format over the months.


### You might find the following interesting
* Functional style python
  * No for loop
  * No while loop
  * Currying, Partial function application
  * Function compisition
* Loading all CSV from subdir in memory, extracting metadata from path/filename using regexp
* Python pivot table through pandas
* Reindexing tables, removing rows/columns or adding empty rows/columns
* Excel production from python
