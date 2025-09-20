#!/usr/bin/env bash
# Download JUnit 5 and required dependencies into ./lib
# Usage: ./download.sh [version]
# Default JUnit version: 5.10.0 (platform 1.10.0)
set -euo pipefail

JUPITER_VERSION="${1:-5.10.0}"
# Map 5.10.x to platform 1.10.x by default
PLATFORM_VERSION="${PLATFORM_VERSION:-1.${JUPITER_VERSION#5.}}"
LIB_DIR="$(cd "$(dirname "$0")/.." && pwd)/lib"

mkdir -p "$LIB_DIR"
cd "$LIB_DIR"

echo "Downloading JUnit Jupiter $JUPITER_VERSION and Platform $PLATFORM_VERSION to $LIB_DIR"

# Core JUnit 5 (Jupiter)
curl -L -o "junit-jupiter-api-${JUPITER_VERSION}.jar"     "https://repo1.maven.org/maven2/org/junit/jupiter/junit-jupiter-api/${JUPITER_VERSION}/junit-jupiter-api-${JUPITER_VERSION}.jar"
curl -L -o "junit-jupiter-engine-${JUPITER_VERSION}.jar"  "https://repo1.maven.org/maven2/org/junit/jupiter/junit-jupiter-engine/${JUPITER_VERSION}/junit-jupiter-engine-${JUPITER_VERSION}.jar"
curl -L -o "junit-jupiter-params-${JUPITER_VERSION}.jar"  "https://repo1.maven.org/maven2/org/junit/jupiter/junit-jupiter-params/${JUPITER_VERSION}/junit-jupiter-params-${JUPITER_VERSION}.jar"

# JUnit Platform
curl -L -o "junit-platform-launcher-${PLATFORM_VERSION}.jar"  "https://repo1.maven.org/maven2/org/junit/platform/junit-platform-launcher/${PLATFORM_VERSION}/junit-platform-launcher-${PLATFORM_VERSION}.jar"
curl -L -o "junit-platform-engine-${PLATFORM_VERSION}.jar"    "https://repo1.maven.org/maven2/org/junit/platform/junit-platform-engine/${PLATFORM_VERSION}/junit-platform-engine-${PLATFORM_VERSION}.jar"
curl -L -o "junit-platform-commons-${PLATFORM_VERSION}.jar"   "https://repo1.maven.org/maven2/org/junit/platform/junit-platform-commons/${PLATFORM_VERSION}/junit-platform-commons-${PLATFORM_VERSION}.jar"

# Additional dependency
curl -L -o "opentest4j-1.3.0.jar" "https://repo1.maven.org/maven2/org/opentest4j/opentest4j/1.3.0/opentest4j-1.3.0.jar"

# API Guardian (transitive dep of junit 5)
curl -L -o "apiguardian-api-1.1.2.jar" "https://repo1.maven.org/maven2/org/apiguardian/apiguardian-api/1.1.2/apiguardian-api-1.1.2.jar"

echo "Done. Files in $LIB_DIR:"
ls -1 *.jar | sed 's/^/  - /'
