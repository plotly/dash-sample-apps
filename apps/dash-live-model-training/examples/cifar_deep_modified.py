# Copyright 2015 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""A deep CIFAR10 classifier using convolutional layers.

Code taken from the official tensorflow guide:
https://www.tensorflow.org/get_started/mnist/pros

The architecture comes from Keras CIFAR10 CNN example:
https://github.com/keras-team/keras/blob/master/examples/cifar10_cnn.py


"""
# Disable linter warnings to maintain consistency with tutorial.
# pylint: disable=invalid-name
# pylint: disable=g-bad-import-order

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys

import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data

# Modified Import
import numpy as np
from sklearn.model_selection import train_test_split
from skimage.transform import rescale
from skimage import color
from tfutils import write_data
from sklearn.preprocessing import OneHotEncoder

FLAGS = None


def deepnn(x):
  """deepnn builds the graph for a deep net for classifying digits.

  Args:
    x: an input tensor with the dimensions (N_examples, 784), where 784 is the
    number of pixels in a standard MNIST image.

  Returns:
    A tuple (y, keep_prob). y is a tensor of shape (N_examples, 10), with values
    equal to the logits of classifying the digit into one of 10 classes (the
    digits 0-9). keep_prob is a scalar placeholder for the probability of
    dropout.
  """
  # Reshape to use within a convolutional neural net.
  # Last dimension is for "features" - there is three here, since images are
  # rgb -- it would be 1 for a grayscale image, 4 for RGBA, etc.
  x_image = tf.reshape(x, [-1, 32, 32, 3])

  # Convolutional layers 1 and 2 - maps 3-color image to 32 feature maps.
  W_conv1 = weight_variable([3, 3, 3, 32])  # 3x3 filters
  b_conv1 = bias_variable([32])
  h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)

  W_conv2 = weight_variable([3, 3, 32, 32])
  b_conv2 = bias_variable([32])
  h_conv2 = tf.nn.relu(conv2d(h_conv1, W_conv2) + b_conv2)

  # Pooling layer - downsamples by 2X.
  h_pool2 = max_pool_2x2(h_conv2)

  # Dropout
  h_pool2_drop = tf.nn.dropout(h_pool2, 0.75)

  # Convolutional layers 3 and 4 - maps 32 feature maps to 64.
  W_conv3 = weight_variable([3, 3, 32, 64])  # 3x3 filters
  b_conv3 = bias_variable([64])
  h_conv3 = tf.nn.relu(conv2d(h_pool2_drop, W_conv3) + b_conv3)

  W_conv4 = weight_variable([3, 3, 64, 64])  # 3x3 filters
  b_conv4 = bias_variable([64])
  h_conv4 = tf.nn.relu(conv2d(h_conv3, W_conv4) + b_conv4)

  # Second pooling layer.
  h_pool4 = max_pool_2x2(h_conv4)

  # Dropout
  h_pool4_drop = tf.nn.dropout(h_pool4, 0.75)

  # Fully connected layer 1 -- after 2 round of downsampling, our 32x32 image
  # is down to 8x8x64 feature maps -- maps this to 512 features.
  W_fc1 = weight_variable([8 * 8 * 64, 512])
  b_fc1 = bias_variable([512])

  h_pool4_flat = tf.reshape(h_pool4_drop, [-1, 8*8*64])
  h_fc1 = tf.nn.relu(tf.matmul(h_pool4_flat, W_fc1) + b_fc1)

  # Dropout - controls the complexity of the model, prevents co-adaptation of
  # features.
  keep_prob = tf.placeholder(tf.float32)
  h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

  # Map the 512 features to 10 classes, one for each digit
  W_fc2 = weight_variable([512, 10])
  b_fc2 = bias_variable([10])

  y_conv = tf.matmul(h_fc1_drop, W_fc2) + b_fc2
  return y_conv, keep_prob


def conv2d(x, W):
  """conv2d returns a 2d convolution layer with full stride."""
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')


def max_pool_2x2(x):
  """max_pool_2x2 downsamples a feature map by 2X."""
  return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                        strides=[1, 2, 2, 1], padding='SAME')


def weight_variable(shape):
  """weight_variable generates a weight variable of a given shape."""
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)


def bias_variable(shape):
  """bias_variable generates a bias variable of a given shape."""
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)


def main(_):
  # Import data
  print("Starting to generate CIFAR10 images.")
  (x_train, y_train), (x_test, y_test) = tf.keras.datasets.cifar10.load_data()
  x_train = np.moveaxis(x_train, 1, 3) / 255.  # Normalize values
  x_train_vec = x_train.reshape(50000, -1)

  x_test = np.moveaxis(x_test, 1, 3) / 255.  # Normalize values
  x_test_vec = x_test.reshape(10000, -1)

  X_train, X_val, y_train, y_val = train_test_split(x_train_vec, y_train, test_size=0.1, random_state=42)
  print("Finished generating CIFAR10 images.")

  # Create the model
  x = tf.placeholder(tf.float32, [None, 32*32*3])

  # Define loss and optimizer
  y_ = tf.placeholder(tf.float32, [None, 10])

  # Build the graph for the deep net
  y_conv, keep_prob = deepnn(x)

  cross_entropy = tf.reduce_mean(
      tf.nn.softmax_cross_entropy_with_logits_v2(labels=y_, logits=y_conv))
  train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)  # RMS is used in keras example, Adam is better
  correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
  accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

  with tf.Session() as sess:
    y_train = OneHotEncoder(sparse=False).fit_transform(y_train)
    y_val = OneHotEncoder(sparse=False).fit_transform(y_val)

    sess.run(tf.global_variables_initializer())
    for i in range(20001):
      start_train = i * 50 % y_train.shape[0]
      end_train = start_train + 50

      start_val = i * 50 % y_val.shape[0]
      end_val = start_val + 50

      batch = (X_train[start_train:end_train], y_train[start_train:end_train])
      batch_val = (X_val[start_val:end_val], y_val[start_val:end_val])

      feed_dict_train = {x: batch[0], y_: batch[1], keep_prob: 1.0}
      feed_dict_val = {x: batch_val[0], y_: batch_val[1], keep_prob: 1.0}
      # Writes data into run log csv file
      write_data(
        accuracy=accuracy,
        cross_entropy=cross_entropy,
        feed_dict_train=feed_dict_train,
        feed_dict_val=feed_dict_val,
        step=i
      )

      if i % 100 == 0:
        train_accuracy = accuracy.eval(feed_dict={
            x: batch[0], y_: batch[1], keep_prob: 1.0})
        print('step %d, training accuracy %g' % (i, train_accuracy))
      train_step.run(feed_dict={x: batch[0], y_: batch[1], keep_prob: 0.5})

    print('test accuracy %g' % accuracy.eval(feed_dict={
        x: x_test_vec, y_: y_test, keep_prob: 1.0}))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--data_dir', type=str,
                      default='/tmp/tensorflow/mnist/input_data',
                      help='Directory for storing input data')
  FLAGS, unparsed = parser.parse_known_args()
  tf.app.run(main=main, argv=[sys.argv[0]] + unparsed)
