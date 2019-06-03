# dash-object-detection

This is a demo of the Dash interactive Python framework developed by [Plotly](https://plot.ly/).

Dash abstracts away all of the technologies and protocols required to build an interactive web-based application and is a simple and effective way to bind a user interface around your Python code. To learn more check out our [documentation](https://plot.ly/dash).

Try out the [demo app here](https://dash-gallery.plotly.host/dash-object-detection/).

![Animated1](images/Screencast.gif)

## Getting Started

### Using the demo

To get started, select a footage you want to view, and choose the display mode (with or without bounding boxes). Then, you can start playing the video, and the visualization will be displayed depending on the current time.

### Running the app locally

First create a virtual environment with conda or venv inside a temp folder, then activate it.

```
virtualenv dash-object-detection

# Windows
dash-object-detection\Scripts\activate
# Or Linux
source venv/bin/activate
```

Clone the git repo, then install the requirements with pip
```
cd dash-object-detection
pip install -r requirements.txt
```

Run the app
```
python app.py
```

## About this app
The videos are displayed using a community-maintained Dash video component. It is made by two Plotly community contributors. You can find the [source code here](https://github.com/SkyRatInd/Video-Engine-Dash).

All videos used are open-sourced under Creative Commons. The [original links can be found here](data/original_footage.md).

### Model
The object detection model is the MobileNet v1, made by Google and trained on the COCO dataset. You can find their implementation on their [official Github repo](https://github.com/tensorflow/models/blob/master/research/slim/nets/mobilenet_v1.md). You are encouraged to try this app with other models.

### Bounding Box Generation
The data displayed in the app are pre-generated for demo purposes. To generate the csv files containing the objects detected for each frame, as well as the output video with bounding boxes, please refer to `utils/generate_video_data.py`. You will need the latest version of tensorflow and OpenCV, as well as the frozen graph `ssd_mobilenet_v1_coco`, that you can [download in the Model Zoo](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/detection_model_zoo.md). Make sure to place the frozen graph inside the same folder as `generate_video_data.py`, i.e. `utils`.

## Built With

* [Dash](https://dash.plot.ly/) - Main server and interactive components
* [Plotly Python](https://plot.ly/python/) - Used to create the interactive plots
* [OpenCV](https://docs.opencv.org/) - Create the video with bounding boxes
* [Tensorflow](https://www.tensorflow.org/api_docs/) - Generate the bounding box data


## Authors

* **Xing Han Lu** - *Initial Work* - [@xhlulu](https://github.com/xhlulu)
* **Yi Cao** - *Restyle* - [@ycaokris](https://github.com/ycaokris)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Screenshots
![Screenshot1](images/Screenshot1.png)

![Screenshot2](images/Screenshot2.png)
