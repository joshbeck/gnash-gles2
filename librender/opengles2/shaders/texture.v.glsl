#version 100
attribute highp vec2 coord2d;
attribute mediump vec2 texcoord;
varying mediump vec2 f_texcoord;

void main(void) {
  gl_Position = vec4(coord2d, 0.0, 1.0);
  f_texcoord = texcoord;
}
