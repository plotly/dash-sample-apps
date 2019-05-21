## Getting Started with the Demo

The purpose of this demo is to show what can be done with the Viewer. Therefore, the trainings shown above are displayed using pre-generated data. To try with your own data, please visit the [project repo](https://github.com/plotly/dash-live-model-training). 

To use the demo, simply choose the model and dataset for which you want to replay the training, __using the two dropdown menus at the top of the page__. For every dataset, we trained a simple 1-layer Neural Network, and a small Convolutional Neural Network that were taken from the official [Tensorflow](https://www.tensorflow.org/tutorials/layers) and [Keras](https://github.com/keras-team/keras/blob/master/examples/cifar10_cnn.py) tutorials.

At the moment, the following datasets have been pre-generated:
* __CIFAR10:__ 50,000 RGB images of size 32x32. It contains 10 common objects. [Link](https://www.cs.toronto.edu/~kriz/cifar.html)
* __MNIST:__ 60,000 grayscale images of size 28x28. It contains handwritten digits from 0 to 9. [Link](http://yann.lecun.com/exdb/mnist/)
* __Fashion MNIST:__ 60,000 grayscale images of size 28x28. It contains 10 commonly found fashion items. [Link](https://github.com/zalandoresearch/fashion-mnist)

## What am I looking at?
The two plots, updated at the interval of your choice, shows two important metrics for training Machine Learning classifiers, i.e. Accuracy and Cross Entropy Loss. Because we feed small mini-batches of 50 examples for every iteration, you will see a lot of variation along the y-axis, which can be made smoother using the sliders. 
 
Accuracy is the fraction of labels correctly classified inside the mini-batch currently used to train/validate the model. In our case, given that each dataset has 10 different labels, an accuracy of 0.1 is equivalent to a random guess, and an accuracy of 1 means the model correctly predicted every single label.

Cross Entropy Loss is the value that you are trying to minimize with your model. It basically indicates how far off our model is from predicting the correct label every time. It is described more in depth in the [Tensorflow Tutorial](https://www.tensorflow.org/versions/r1.0/get_started/mnist/beginners#training).

## What does the app do?
For the majority of Deep Learning models, it is extremely helpful to keep track of the accuracy and loss as it is training. At the moment, the best application to do that is [Tensorboard](https://www.tensorflow.org/programmers_guide/summaries_and_tensorboard), which is a collection of visualization tools (metrics plots, image examples, graph representation, weight histogram, etc.) useful to debug and monitor the training of your model.

_Dash's Live Model Training Viewer_ is a compact visualization app that monitors core metrics of your __Tensorflow model__ during training. It complements the Tensorboard by offering the following:
* __Real-time visualization__: The app is designed to visualize your metrics as they are updated inside your model.
* __Small and Lightweight__: The viewer loads a small number of important visualization, so that it loads and runs quickly.
* __Simple to use__: For simpler tensorflow models, all you need to do is to call `add_eval` to add the accuracy and cross entropy operations in the graph, and generate a log of the metrics using `write_data`. Both functions are inside `tfutils.py`, and examples are included in the `examples` directory.
* __Easy to modify__: The app is stored inside one module, and is written in under 400 lines. You can quickly modify and improve the app without breaking anything.
* __Plotly Graphs and Dash Integration__: Easily integrate the app into more complex Dash Apps, and includes all the tools found in Plotly graphs.

Here is a flowchart of the process from model training to visualization:
![flowchart](https://raw.githubusercontent.com/plotly/dash-live-model-training/master/images/flowchart.png)

At the moment, the logging only works for iterative Tensorflow models. We are planning to extend it for PyTorch. You are encouraged to port the logging function (which is a simple csv logging) to Keras, Tensorflow's high-level API, MXNet, etc.

You can find the code for the app, and examples of using the app in the [repository](https://github.com/plotly/dash-live-model-training).