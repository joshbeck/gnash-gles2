#version 100
attribute highp vec2 coord2d;
attribute lowp vec4 v_color;
varying lowp vec4 f_color;

void main(void) {
  gl_Position = vec4(coord2d, 0.0, 1.0);
  f_color = v_color;
}
