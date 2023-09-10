using Plots
using LaTeXStrings

x = 0:0.01:10
y = x .^ 2
plot_variable = plot(x, y, label=L"y = x^2")
display(plot_variable)
print("Hello world!!")
 