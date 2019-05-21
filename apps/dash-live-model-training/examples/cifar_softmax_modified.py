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
"""A very simple CIFAR10 classifier.

This code was modified from the MNIST beginner tutorial found here:
https://www.tensorflow.org/get_started/mnist/beginners

Accuracy on test set is 26.74%. Simple ConvNet models can achieve over 70% accuracy:
https://github.com/keras-team/keras/blob/master/examples/cifar10_cnn.py

Whether the low accuracy is caused by an error or the simplicity of the classifier is unknown. It is encouraged to
report errors within this code.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys

import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data

# Custom Imports
import numpy as np
from sklearn.model_selection import train_test_split
from skimage.transform import rescale
from skimage import color
from tfutils import add_eval, write_data

FLAGS = None


def main(_):
  # Import data
  print("Starting to generate CIFAR10 images.")
  (x_train, y_train), (x_test, y_test) = tf.keras.datasets.cifar10.load_data()
  x_train = np.moveaxis(x_train, 1, 3) / 255.  # Normalize values
  x_train_vec = x_train.reshape(50000, -1)

  y_train = np.squeeze(y_train)
  y_test = np.squeeze(y_test)

  x_test = np.moveaxis(x_test, 1, 3) / 255.  # Normalize values
  x_test_vec = x_test.reshape(10000, -1)

  X_train, X_val, y_train, y_val = train_test_split(x_train_vec, y_train, test_size=0.1, random_state=42)
  print("Finished generating CIFAR10 images.")

  # Create the model
  x = tf.placeholder(tf.float32, [None, 3*32*32])
  W = tf.Variable(tf.zeros([3*32*32, 10]))
  b = tf.Variable(tf.zeros([10]))
  y = tf.matmul(x, W) + b

  # Define loss and optimizer
  y_ = tf.placeholder(tf.int64, [None])

  # The raw formulation of cross-entropy,
  #
  #   tf.reduce_mean(-tf.reduce_sum(y_ * tf.log(tf.nn.softmax(y)),
  #                                 reduction_indices=[1]))
  #
  # can be numerically unstable.
  #
  # So here we use tf.losses.sparse_softmax_cross_entropy on the raw
  # outputs of 'y', and then average across the batch.
  cross_entropy = tf.losses.sparse_softmax_cross_entropy(labels=y_, logits=y)
  train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)

  # Add accuracy and cross entropy to the graph using util function
  accuracy, cross_entropy = add_eval(y, y_)

  sess = tf.InteractiveSession()
  tf.global_variables_initializer().run()
  # Train
  for i in range(20001):
    start_train = i * 100 % y_train.shape[0]
    end_train = start_train + 100

    start_val = i * 100 % y_val.shape[0]
    end_val = start_val + 100

    batch = (X_train[start_train:end_train], y_train[start_train:end_train])
    batch_val = (X_val[start_val:end_val], y_val[start_val:end_val])

    feed_dict_train = {x: batch[0], y_: batch[1]}
    feed_dict_val = {x: batch_val[0], y_: batch_val[1]}
    # Writes data into run log csv file
    write_data(
        accuracy=accuracy,
        cross_entropy=cross_entropy,
        feed_dict_train=feed_dict_train,
        feed_dict_val=feed_dict_val,
        step=i
    )
    sess.run(train_step, feed_dict={x: batch[0], y_: batch[1]})

  # Test trained model
  correct_prediction = tf.equal(tf.argmax(y, 1), y_)
  accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
  print(sess.run(
      accuracy, feed_dict={
          x: x_test_vec,
          y_: y_test
      }))


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument(
      '--data_dir',
      type=str,
      default='/tmp/tensorflow/mnist/input_data',
      help='Directory for storing input data')
  FLAGS, unparsed = parser.parse_known_args()
  tf.app.run(main=main, argv=[sys.argv[0]] + unparsed)
