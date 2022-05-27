from datetime import datetime as dt
import numpy as np

from constants import totalList

def get_selection(month, day, selection):
    """
        Get the amount of rides per hour based on the time selected
        This also higlights the color of the histogram bars based on
        if the hours are selected
    """
    xVal, yVal, xSelected = [], [], []
    colorVal = [
        "#F4EC15",
        "#DAF017",
        "#BBEC19",
        "#9DE81B",
        "#80E41D",
        "#66E01F",
        "#4CDC20",
        "#34D822",
        "#24D249",
        "#25D042",
        "#26CC58",
        "#28C86D",
        "#29C481",
        "#2AC093",
        "#2BBCA4",
        "#2BB5B8",
        "#2C99B4",
        "#2D7EB0",
        "#2D65AC",
        "#2E4EA4",
        "#2E38A4",
        "#3B2FA0",
        "#4E2F9C",
        "#603099",
    ]

    # Put selected times into a list of numbers xSelected
    xSelected.extend([int(x) for x in selection])

    for i in range(24):
        # If bar is selected then color it white
        if i in xSelected and len(xSelected) < 24:
            colorVal[i] = "#FFFFFF"
        xVal.append(i)
        # Get the number of rides at a particular time
        yVal.append(len(totalList[month][day][totalList[month][day].index.hour == i]))
    return [np.array(xVal), np.array(yVal), np.array(colorVal)]


def getLatLonColor(selectedData, month, day):
    " Get the Coordinates of the chosen months, dates and times "
    listCoords = totalList[month][day]

    # No times selected, output all times for chosen month and date
    if selectedData is None or len(selectedData) is 0:
        return listCoords
    listStr = "listCoords["
    for time in selectedData:
        if selectedData.index(time) is not len(selectedData) - 1:
            listStr += "(totalList[month][day].index.hour==" + str(int(time)) + ") | "
        else:
            listStr += "(totalList[month][day].index.hour==" + str(int(time)) + ")]"
    return eval(listStr)


def total_rides_calculation(date_picked, bars_selected):
    
    firstOutput = ""

    if bars_selected is not None or len(bars_selected) != 0:
        date_temp = dt.strptime(date_picked, "%Y-%m-%d")
        totalInSelection = 0
        for x in bars_selected:
            totalInSelection += len(
                totalList[date_temp.month - 4][date_temp.day - 1][
                    totalList[date_temp.month - 4][date_temp.day - 1].index.hour
                    == int(x)
                ]
            )
        firstOutput = "Total rides in selection: {:,d}".format(totalInSelection)

    if (bars_selected is None
        or bars_selected is None
        or len(bars_selected)==24
        or len(bars_selected)==0
    ):
        return firstOutput, (date_picked, " - showing hour(s): All")

    holder = sorted([int(x) for x in bars_selected])
    if holder == list(range(min(holder), max(holder) + 1)):
        return (
            firstOutput,
            (
                date_picked,
                " - showing hour(s): ",
                holder[0],
                "-",
                holder[len(holder) - 1],
            ),
        )
    holder_to_string = ", ".join(str(x) for x in holder)
    return firstOutput, (date_picked, " - showing hour(s): ", holder_to_string)