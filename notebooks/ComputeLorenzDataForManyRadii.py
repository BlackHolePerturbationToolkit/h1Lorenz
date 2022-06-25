import subprocess

r0i = 9
while r0i <= 30:
	r0f = 9
	while r0f <= 9:
		subprocess.run(["math", "-script", "Circular_Lorenz_output_horiz_and_inf_fields_script.m", str(r0i), str(r0f)])
		r0f += 1
	r0i += 1
