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
"""A very simple MNIST classifier.

See extensive documentation at
https://www.tensorflow.org/get_started/mnist/beginners
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys

import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data

from tfutils import add_eval, write_data

FLAGS = None
DATA = "MNIST"


def main(_):
    # Import data
    if DATA == "MNIST":
        mnist = input_data.read_data_sets(FLAGS.data_dir)
    elif DATA == "FASHION":
        mnist = input_data.read_data_sets(
            "data/fashion",
            source_url="http://fashion-mnist.s3-website.eu-central-1.amazonaws.com/",
        )
    # Create the model
    x = tf.placeholder(tf.float32, [None, 784])
    W = tf.Variable(tf.zeros([784, 10]))
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

    ################################## MODIFIED CODE BELOW ##################################
    accuracy, cross_entropy = add_eval(y, y_)
    ################################## MODIFIED CODE ABOVE ##################################

    sess = tf.InteractiveSession()
    tf.global_variables_initializer().run()
    # Train
    for i in range(10001):
        batch_xs, batch_ys = mnist.train.next_batch(100)

        ################################## MODIFIED CODE BELOW ##################################
        batch = mnist.train.next_batch(100)
        batch_val = mnist.validation.next_batch(100)
        feed_dict_train = {x: batch[0], y_: batch[1]}
        feed_dict_val = {x: batch_val[0], y_: batch_val[1]}
        # Writes data into run log csv file
        write_data(
            accuracy=accuracy,
            cross_entropy=cross_entropy,
            feed_dict_train=feed_dict_train,
            feed_dict_val=feed_dict_val,
            step=i,
        )
        ################################## MODIFIED CODE ABOVE ##################################

        sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})

    # Test trained model
    correct_prediction = tf.equal(tf.argmax(y, 1), y_)
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    print(sess.run(accuracy, feed_dict={x: mnist.test.images, y_: mnist.test.labels}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data_dir",
        type=str,
        default="/tmp/tensorflow/mnist/input_data",
        help="Directory for storing input data",
    )
    FLAGS, unparsed = parser.parse_known_args()
    tf.app.run(main=main, argv=[sys.argv[0]] + unparsed)
