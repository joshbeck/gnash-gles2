#version 100
varying lowp vec4 f_color;

void main(void) {
    gl_FragColor = vec4(f_color.r, f_color.g, f_color.b, f_color.a);
}
