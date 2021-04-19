// Sass configuration
const gulp = require('gulp');
const sass = require('gulp-sass');

// compile scss into css
function style() {
  // outputStyle props: nested, expanded, compact, compressed
  return gulp.src('../scss/**/*.scss')
    .pipe(sass({outputStyle: 'compressed'}).on('error', sass.logError))
    .pipe(gulp.dest('../assets/css'));
};

function watch(){
  gulp.watch('../scss/**/*.scss', style);
};

exports.style = style;
exports.watch = watch;