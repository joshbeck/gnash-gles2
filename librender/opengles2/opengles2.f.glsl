#version 100
varying mediump vec2 f_texcoord;
varying lowp vec4 f_color;
uniform sampler2D texture;

void main(void) {
  if(!texture) {
    gl_FragColor = vec4(f_color.r, f_color.g, f_color.b, f_color.a);
  } else {
    gl_FragColor = texture2D(texture, f_texcoord);
  }
}
