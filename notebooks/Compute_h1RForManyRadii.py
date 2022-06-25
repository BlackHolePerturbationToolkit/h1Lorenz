import subprocess

r0s = [9.1, 9.2]
for r0 in r0s:
	subprocess.run(["math", "-script", "Calc_h1R_ab_for_script.m", str(r0)])
