#pragma once

#pragma once

#include <cstddef>
#include <filesystem>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace HighFive {

class DataSpace {
public:
  DataSpace() : m_size(0) {}
  DataSpace(std::initializer_list<std::size_t> dims) : m_size(product(dims)) {}
  explicit DataSpace(const std::vector<std::size_t> &dims)
      : m_size(product_vector(dims)) {}

  template <typename T> static DataSpace From(const T &value) {
    if constexpr (requires { value.size(); }) {
      return DataSpace({static_cast<std::size_t>(value.size())});
    } else {
      return DataSpace({1});
    }
  }

  [[nodiscard]] std::size_t size() const { return m_size; }

private:
  static std::size_t product(std::initializer_list<std::size_t> dims) {
    std::size_t result = 1;
    if (dims.size() == 0) {
      return 0;
    }
    for (auto dim : dims) {
      result *= dim;
    }
    return result;
  }

  static std::size_t product_vector(const std::vector<std::size_t> &dims) {
    if (dims.empty()) {
      return 0;
    }
    std::size_t result = 1;
    for (auto dim : dims) {
      result *= dim;
    }
    return result;
  }

  std::size_t m_size;
};

namespace detail {

struct DataBufferBase {
  virtual ~DataBufferBase() = default;
};

template <typename T> struct DataBuffer : DataBufferBase {
  std::vector<T> data;
};

struct FileStorage {
  std::unordered_map<std::string, std::shared_ptr<DataBufferBase>> datasets;
};

inline auto &file_registry() {
  static std::unordered_map<std::string, std::shared_ptr<FileStorage>> registry;
  return registry;
}

inline std::shared_ptr<FileStorage> get_storage(const std::string &filename,
                                                bool overwrite) {
  auto &registry = file_registry();
  if (overwrite) {
    auto storage = std::make_shared<FileStorage>();
    registry[filename] = storage;
    return storage;
  }

  auto it = registry.find(filename);
  if (it == registry.end() || it->second == nullptr) {
    throw std::runtime_error("HighFive stub: file not found: " + filename);
  }
  return it->second;
}

inline std::string join_paths(const std::string &base,
                              const std::string &child) {
  if (base.empty()) {
    return child;
  }
  if (base.back() == '/') {
    return base + child;
  }
  return base + "/" + child;
}

} // namespace detail

class DataSet {
public:
  DataSet() = default;

  template <typename T> void write(const std::vector<T> &data) {
    if (!m_storage) {
      throw std::runtime_error("HighFive stub: dataset without storage");
    }
    auto buffer = std::make_shared<detail::DataBuffer<T>>();
    buffer->data = data;
    m_storage->datasets[m_path] = buffer;
  }

  template <typename T> void read(std::vector<T> &out) const {
    if (!m_storage) {
      throw std::runtime_error("HighFive stub: dataset without storage");
    }
    auto it = m_storage->datasets.find(m_path);
    if (it == m_storage->datasets.end()) {
      throw std::runtime_error("HighFive stub: dataset not found: " + m_path);
    }
    auto buffer = std::dynamic_pointer_cast<detail::DataBuffer<T>>(it->second);
    if (!buffer) {
      throw std::runtime_error("HighFive stub: dataset type mismatch");
    }
    out = buffer->data;
  }

private:
  friend class File;
  friend class Group;

  DataSet(std::shared_ptr<detail::FileStorage> storage, std::string path)
      : m_storage(std::move(storage)), m_path(std::move(path)) {}

  std::shared_ptr<detail::FileStorage> m_storage;
  std::string m_path;
};

template <typename T> class Attribute {
public:
  Attribute() = default;
  template <typename Value> void write(const Value &) const {}
};

class Group {
public:
  Group() = default;

  Group createGroup(const std::string &name) const {
    return Group(m_storage, detail::join_paths(m_path, name));
  }

  template <typename T>
  DataSet createDataSet(const std::string &name,
                        const std::vector<T> &data) const {
    auto dataset = DataSet(m_storage, detail::join_paths(m_path, name));
    dataset.write(data);
    return dataset;
  }

  template <typename T>
  DataSet createDataSet(const std::string &name, const DataSpace &space) const {
    std::vector<T> data(space.size());
    return createDataSet(name, data);
  }

  template <typename T>
  Attribute<T> createAttribute(const std::string &, const DataSpace &) const {
    return Attribute<T>();
  }

private:
  friend class File;

  Group(std::shared_ptr<detail::FileStorage> storage, std::string path)
      : m_storage(std::move(storage)), m_path(std::move(path)) {}

  std::shared_ptr<detail::FileStorage> m_storage;
  std::string m_path;
};

class File {
public:
  struct OverwriteTag {};
  struct ReadOnlyTag {};

  static constexpr OverwriteTag Overwrite{};
  static constexpr ReadOnlyTag ReadOnly{};

  template <typename Path>
  File(Path path, OverwriteTag) : m_filename(path_to_string(path)) {
    m_storage = detail::get_storage(m_filename, true);
  }

  template <typename Path>
  File(Path path, ReadOnlyTag) : m_filename(path_to_string(path)) {
    m_storage = detail::get_storage(m_filename, false);
  }

  Group createGroup(const std::string &name) const {
    return Group(m_storage, name);
  }

  DataSet getDataSet(const std::string &path) const {
    return DataSet(m_storage, path);
  }

private:
  template <typename Path> static std::string path_to_string(const Path &path) {
    if constexpr (std::is_same_v<std::decay_t<Path>, std::filesystem::path>) {
      return path.string();
    } else {
      return std::string(path);
    }
  }

  std::shared_ptr<detail::FileStorage> m_storage;
  std::string m_filename;
};

} // namespace HighFive
