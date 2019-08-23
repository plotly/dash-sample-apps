ternary_color = list('#FFDF7F', '#FFDF7F', '#FFDF7F',
                 '#FFF9DC', '#FFF9DC', '#FFF9DC', '#FFF9DC', '#FFF9DC', '#FFF9DC',
                 '#FFF2AD','#FFF2AD', '#FFF2AD',
                 '#FDFFE9', '#FDFFE9', '#FDFFE9', '#FDFFE9')

# Colors for legend
colors = list(
  "#001f3f",
  "#0074d9",
  "#3d9970",
  "#111111",
  "#01ff70",
  "#ffdc00",
  "#ff851B",
  "#ff4136",
  "#85144b",
  "#f012be",
  "#b10dc9",
  "#AAAAAA",
  "#111111"
)

# Ternary axis contour
ternary_contour = list(
  'Carbonate Dominated Lith'= list(
    list('Quartz'= 0, 'Carbonate'= 100, 'Clay'= 0),
    list('Quartz'= 20, 'Carbonate'= 80, 'Clay'= 0),
    list('Quartz'= 0, 'Carbonate'= 80, 'Clay'= 20)
    ),
  'Silica Dominated Lith'= list(
    list('Quartz'= 100, 'Carbonate'= 0, 'Clay'= 0),
    list('Quartz'= 80, 'Carbonate'= 20, 'Clay'= 0),
    list('Quartz'= 80, 'Carbonate'= 0, 'Clay'= 20)
    ),
  'Clay Dominated Lith'= list(
    list('Quartz'= 0, 'Carbonate'= 20, 'Clay'= 80),
    list('Quartz'= 0, 'Carbonate'= 0, 'Clay'= 100),
    list('Quartz'= 20, 'Carbonate'= 0, 'Clay'= 80)
    ),
  'Silica rich Carbonate Mudstone'= list(
    list('Quartz'= 20, 'Carbonate'= 80, 'Clay'= 0),
    list('Quartz'= 10, 'Carbonate'= 80, 'Clay'= 10),
    list('Quartz'= 40, 'Carbonate'= 50, 'Clay'= 10),
    list('Quartz'= 50, 'Carbonate'= 50, 'Clay'= 0)
    ),
  'Silica rich Argillaceous Mudstone'= list(
    list('Quartz'= 40, 'Carbonate'= 10, 'Clay'= 50),
    list('Quartz'= 10, 'Carbonate'= 10, 'Clay'= 80),
    list('Quartz'= 20, 'Carbonate'= 0, 'Clay'= 80),
    list('Quartz'= 50, 'Carbonate'= 0, 'Clay'= 50)
    ),
  'Carbonate-rich Argillaceous Mudstone'= list(
    list('Quartz'= 0, 'Carbonate'= 50, 'Clay'= 50),
    list('Quartz'= 0, 'Carbonate'= 20, 'Clay'= 80),
    list('Quartz'= 10, 'Carbonate'= 10, 'Clay'= 80),
    list('Quartz'= 10, 'Carbonate'= 40, 'Clay'= 50)
    ),
  'Carbonate rich Silliceous Mudstone'= list(
    list('Quartz'= 80, 'Carbonate'= 10, 'Clay'= 10),
    list('Quartz'= 50, 'Carbonate'= 10, 'Clay'= 40),
    list('Quartz'= 50, 'Carbonate'= 0, 'Clay'= 50),
    list('Quartz'= 80, 'Carbonate'= 0, 'Clay'= 20)
    ),
  'Clay rich Silliceous Mudstone'= list(
    list('Quartz'= 80, 'Carbonate'= 20, 'Clay'= 0),
    list('Quartz'= 80, 'Carbonate'= 10, 'Clay'= 10),
    list('Quartz'= 50, 'Carbonate'= 40, 'Clay'= 10),
    list('Quartz'= 50, 'Carbonate'= 50, 'Clay'= 0)
    ),
  'Clay-rich Carbonate Mudstone'= list(
    list('Quartz'= 10, 'Carbonate'= 80, 'Clay'= 10),
    list('Quartz'= 0, 'Carbonate'= 80, 'Clay'= 20),
    list('Quartz'= 0, 'Carbonate'= 50, 'Clay'= 50),
    list('Quartz'= 10, 'Carbonate'= 50, 'Clay'= 40)
    ),
  'Carbonate Silliceous Mudstone'= list(
    list('Quartz'= 50, 'Carbonate'= 50, 'Clay'= 0),
    list('Quartz'= 30, 'Carbonate'= 50, 'Clay'= 20),
    list('Quartz'= 50, 'Carbonate'= 30, 'Clay'= 20)
    ),
  'Argillaceous Carbonate Mudstone'= list(
    list('Quartz'= 20, 'Carbonate'= 50, 'Clay'= 30),
    list('Quartz'= 0, 'Carbonate'= 50, 'Clay'= 50),
    list('Quartz'= 20, 'Carbonate'= 30, 'Clay'= 50)
    ),
  'Argillaceous Sillicaous Mudstone'= list(
    list('Quartz'= 50, 'Carbonate'= 20, 'Clay'= 30),
    list('Quartz'= 30, 'Carbonate'= 20, 'Clay'= 50),
    list('Quartz'= 50, 'Carbonate'= 0, 'Clay'= 50)
    ),
  'Mixed Silliceous Mudstone'= list(
    list('Quartz'= 80, 'Carbonate'= 10, 'Clay'= 10),
    list('Quartz'= 50, 'Carbonate'= 40, 'Clay'= 10),
    list('Quartz'= 50, 'Carbonate'= 10, 'Clay'= 40)
    ),
  'Mixed Mudstone'= list(
    list('Quartz'= 30, 'Carbonate'= 50, 'Clay'= 20),
    list('Quartz'= 20, 'Carbonate'= 50, 'Clay'= 30),
    list('Quartz'= 20, 'Carbonate'= 30, 'Clay'= 50),
    list('Quartz'= 30, 'Carbonate'= 20, 'Clay'= 50),
    list('Quartz'= 50, 'Carbonate'= 20, 'Clay'= 30),
    list('Quartz'= 50, 'Carbonate'= 30, 'Clay'= 20)
    ),
  'Mixed Argillaceous Mudstone'= list(
    list('Quartz'= 40, 'Carbonate'= 10, 'Clay'= 50),
    list('Quartz'= 10, 'Carbonate'= 40, 'Clay'= 50),
    list('Quartz'= 10, 'Carbonate'= 10, 'Clay'= 80)
    ),
  'Mixed Carbonate Mudstone'= list(
    list('Quartz'= 40, 'Carbonate'= 50, 'Clay'= 10),
    list('Quartz'= 10, 'Carbonate'= 80, 'Clay'= 10),
    list('Quartz'= 10, 'Carbonate'= 50, 'Clay'= 40)
    )
)