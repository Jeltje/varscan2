# Definitions
build_path = build/
runtime_path = runtime/
build_output = ${runtime_path}/exec/varscan.jar
runtime_fullpath = $(realpath ${runtime_path})
build_tool = runtime-container.DONE
nametag = jeltje/varscan

# Steps
all: ${build_output} ${build_tool}

${build_output}: ${build_path}/Dockerfile
	cd ${build_path} && docker build -t varscanbuild .
	docker run -v ${runtime_fullpath}:/data varscanbuild cp -rp exec /data

${build_tool}: ${build_output} ${runtime_path}/Dockerfile
	cd ${runtime_path} && docker build -t ${nametag} .
	docker rmi -f varscanbuild
	touch ${build_tool}

#push: all
#	# Requires ~/.dockercfg
#	docker push ${nametag}
