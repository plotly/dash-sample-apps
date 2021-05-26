# Fashion MNIST Explorer with T-SNE and CNN

This app displays the output of two models on the [fashion MNIST dataset](https://github.com/zalandoresearch/fashion-mnist), a collection of 70,000 low resolution images (28 pixel by 28 pixel) belonging to 10 different classes of clothing items.

### Model 1: T-SNE

T-SNE is a dimensionality reduction algorithm that allows to embed high dimensionality data (such as images) onto few dimensions as a means of visualizing the similarity between different data points. The model is pre-computed, and the embeddings were generated using [RAPIDS](https://rapids.ai/about.html). Using a single, modest GPU (GeForce GTX 1060), the embeddings take about 3 seconds to generate for the 70,000, unreduced fashion MNIST images. The code to generate the embeddings can be found in `tsne.py`.

### Model 2: CNN
The fashion MNIST images were split into a training and testing set, and fed into a small convolutional neural network (cnn). This cnn was fit using keras, and training over 5 epochs took under two minutes with a GPU, achieving about 89% accuracy on the testing data. The code to fit the model can be found in `cnn.py`.

## What this app does and how to use it:

Users can choose to display the TSNE embeddings of the training images, testing images, or both over two dimensions. The color of the markers in the scatter correspond to the true class that each image belongs to. By hovering over points, users see a thumbnail of the image represented by that point. Clicking on a point feeds the image it represents into the cnn, and outputs a prediction. The prediction is displayed along with a percentage of certainty.

Users can also upload their own images, which get transformed and classified by the cnn. Some sample images are provided in the repo.

