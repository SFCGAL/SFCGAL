# SFCGAL Continuous Integration Status

## Overview of CI Platforms

| Platform | Badge | Status | Environment |
|:---------|:------|:-------|:------------|
| **GitLab** | [![GitLab pipeline status](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | 🚀 Main Pipeline | All environments |

## Detailed CI Status

### 🦊 GitLab Pipelines

#### 📊 Sonar Analysis
| Job | Badge | Environment |
|:----|:------|:------------|
| Sonar | [![Sonar](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg?job=sonar-build-test)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | Debian |

#### 🧪 Platform Testing
| Job | Badge | Environment |
|:----|:------|:------------|
| Platform Tests | [![Tests](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg?job=platform-build-test)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | Multiple Environments |

##### 🐧 Debian Environments
| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| Debian Stable | Clang | CGAL 6.1.3 |
| Debian Testing | Clang | CGAL 6.1.3 |
| Debian Stable | Clang | CGAL 6.2.1 |
| Debian Testing | Clang | CGAL 6.2.1 |
| Debian Stable | GCC | CGAL 6.1.3 |
| Debian Testing | GCC | CGAL 6.1.3 |
| Debian Stable | GCC | CGAL 6.2.1 |
| Debian Testing | GCC | CGAL 6.2.1 |

##### 🎩 Fedora Environments
| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| Fedora 43 | GCC | CGAL 6.1.3 |
| Fedora 43 | GCC | CGAL 6.2.1 |
| Fedora 43 | Clang | CGAL 6.1.3 |
| Fedora 43 | Clang | CGAL 6.2.1 |
| Fedora 44 | GCC | CGAL 6.1.3 |
| Fedora 44 | GCC | CGAL 6.2.1 |
| Fedora 44 | Clang | CGAL 6.1.3 |
| Fedora 44 | Clang | CGAL 6.2.1 |

##### 🦎 OpenSUSE Environment
| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| OpenSUSE Tumbleweed | GCC | CGAL System* |

#### 🐳 Docker Builds
| Job | Badge | Environment |
|:----|:------|:------------|
| Docker | [![Docker](https://gitlab.com/sfcgal/SFCGAL/badges/master/pipeline.svg?job=docker)](https://gitlab.com/sfcgal/SFCGAL/-/commits/master) | Multiple Environments |

| Environment | Type |
|:------------|:-----|
| Docker Debian | 🐧 Linux |
| Docker Windows | 🪟 Windows |

### 😺 GitHub CI

#### 🔄 GitHub Actions

##### 🔍 CodeQL Analysis
| Platform | Badge | Environment |
|:---------|:------|:------------|
| CodeQL | [![CodeQL](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/codeql.yml) | Debian GCC CGAL 6.2.1 |

##### 📦 VCPKG Build
| Platform | Badge | Environment |
|:---------|:------|:------------|
| VCPKG | [![Build with vcpkg](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/vcpkg.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/vcpkg.yml) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| 🐧 Ubuntu Latest | System | CGAL System* |
| 🍎 macOS Latest | System | CGAL System* |
| 🪟 Windows Latest | System | CGAL System* |

##### 🪟  MSYS2 Build
| Platform | Badge | Environment |
|:---------|:------|:------------|
| MSYS2 | [![MSYS2](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/msys.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/msys.yml) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| mingw64 | System | CGAL System* |
| ucrt64 | System | CGAL System* |
| clang64 | System | CGAL System* |

##### 🌐 Cross Platform Actions
| Platform | Badge | Environment |
|:---------|:------|:------------|
| Cross Platform | [![CI](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/bsd.yml/badge.svg?branch=master)](https://github.com/Oslandia/SFCGAL_CI/actions/workflows/bsd.yml) | Multiple Environments |

| Environment | Compiler | CGAL Version |
|:------------|:---------|:-------------|
| 😈 FreeBSD 14.5 (quarterly) | Clang | CGAL 6.2.1 |
| 🚩 NetBSD 11.0 | GCC | CGAL 6.2.1 |
| 🐡 OpenBSD 7.9 | Clang | CGAL 6.2.1 |
| 🍎 macOS | AppleClang | CGAL System* |

\* *CGAL System refers to the CGAL version available through the system's package manager*
