# Animal Study Browser


> This is the implementation of the Animal Study Browser app for Dash for R.  
The original Dash for Python version is located [here](https://github.com/shday/studybrowser).

This app displays the results of a study comparing several treatment 
groups and optionally calculates p-values. Data from one or more
studies are loaded from a csv file containing the following 
column headers (case sensitive):

* *study_id*
* *group_id*
* *group_type*
* *reading_value*

The file should have one row for each subject. If a group has
group_type "control", it will be compared to the other groups using a t-test.
The group_type "reference" will suppress this calculation. 

The app will also recognize the following columns, providing
some additional functionality:

* *subject_id* - displayed on hover-over of data points
* *test_article* - displayed in the study selection drop-down
* *group_name* - replaces group_id on the x-axis
* *reading_name* - y-axis title 


### Animal studies in drug discovery
Animal studies are often used to determine if a drug candidate (i.e., test article) 
has the desired effect in disease model. For example, a study 
could measure the ability of a compound to reduce tumor size in a mouse cancer 
model or control blood sugar in diabetic rats. Studies are almost always designed to
compare against untreated controls and test for statistical significance. 




## Running the app locally



Clone the git repo, then install the requirements 
```
git clone https://github.com/plotly/dash-sample-apps/dashr-studybrowser.git
cd dash-sample-apps/dashr-studybrowser

```
Download the script and run init.R

```
Rscript init.R
```
Run the app
```
Rscript app.R
```

Open a browser at http://127.0.0.1:8050

## Screenshots


