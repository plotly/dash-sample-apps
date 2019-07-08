### How to use the demo app

To get started, choose a dataset you want to visualize. Here are the image datasets available:

- **MNIST Digits:** Handwritten digits between 0 and 9, of shape 28x28.
- **Fashion MNIST:** Ten common pieces of clothing, of shape 28x28.
- **CIFAR 10:** Ten common animals and vehicles, of shape 32x32.

When the scatter plot appears on the graph, you can see the original image by clicking on a data point.

Alternatively, you can explore the [GloVe Word Vectors datasets](https://nlp.stanford.edu/projects/glove/), which are encoded vectors of large collection of texts from Wikipedia, Twitter, and acquired through Web Crawlers. You will be offered two ways of displaying the word embeddings:

- **Regular:** All 3000 samples will be shown in scatter dots.
- **Nearest Neighbors:** The 100 closest neighbors of the word selected, in term of Euclidean distances. To select a word, using the dropdown menu for this purpose. The neighbors will be shown in text.

Upon clicking a data point, you will be able to see the 5 closest neighbors of the word you clicked.

### How does t-SNE work?

Images can be seen as long vectors of hundred or thousands of dimensions, each dimension representing one shade of color, or of gray (if the image is black & white). For example, a 28x28 image (such as MNIST) can be unrolled into one 784-dimensional vector. Kept this way, it would be extremely hard to visualize our dataset, especially if it contains tens of thousands of samples; in fact, the only way would be to go through each and every image, and keep track of how they are all written.

The t-SNE algorithm solves this problem by reducing the number of dimensions of your datasets, so that you can visualize it in low-dimensional space, i.e. in 2D or 3D. For each data point, you will now have a position on your 3D plot, which can be compared with other data points to understand how close or far apart they are from each other.

For the MNIST Digits dataset, we can see that sometimes there are a few digits clustered together, with a few outliers probably caused by poor handwriting. When that happens, click on the images to see what kind of problem it might have; if none, maybe try to change the input parameters to find a tighter clustering.

### Choosing the right parameters

The quality of a t-SNE visualization depends heavily on the input parameters when you train the algorithm. Each parameter has a great impact on how well each group of data will be clustered. Here is what you should know for each of them:

- **Number of Iterations:** This is how many steps you want to run the algorithm. A higher number of iterations often gives better visualizations, but more time to train.
- **Perplexity:** This is a value that influences the number of neighbors that are taken into account during the training. According to the [original paper](https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf), the value should be between 5 and 50.
- **Learning Rate:** This value determines how much weight we should give to the updates given by the algorithm at each step. It is typically between 10 and 1000.
- **Initial PCA Dimensions:** Because the number of dimensions of the original data might be very big, we use another [dimensionality reduction technique called PCA](https://en.wikipedia.org/wiki/Principal_component_analysis) to first reduce the dataset to a smaller space, and then apply t-SNE on that space. Initial PCA dimensions of 50 has shown good experimental results in the original paper.

### GloVe & Word Embeddings

Image datasets are straightforward, since the high dimensionality comes from the numbers of pixels in the image. However, what if we want to visualize how close or far two words are, with respect to how often to appear with each other? One way to do so is to visualize word embeddings, which are large collections of text (e.g. Wikipedia Articles, User Tweets, Web Scraped texts) reduced to smaller dimensions (50-d, 100-d, 200-d, etc.) that contain meaningful and compact information about how often words appear with each other, with respect to how often they are used. We used [GloVe embeddings, which are developed by the Stanford NLP group](). Andrew Ng offers intuitive and insightful explanations about Word embeddings through his free [Sequence Models Course](https://www.coursera.org/learn/nlp-sequence-models/lecture/qHMK5/using-word-embeddings).
