/*----------------------------------------------------------------------------*/
#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include <stdexcept>

namespace gecko {
    // 1. Forward declaration of the friend class
    // 1. Forward declaration du template
    template<typename TopologyTraits>
    class UnstructuredMesh;

    /**
     * @brief Type-erased interface implemented by every Variable<T>, allowing
     * VariableRegistry to resize and manipulate variables without knowing their type.
     */
    /*----------------------------------------------------------------------------*/
    struct IVariable {
        virtual ~IVariable() = default;

        /** @brief Gets the number of items currently stored. @return The number of items. */
        [[nodiscard]] virtual size_t size() const = 0;

        /**
         * @brief Resizes the underlying storage.
         * @param new_size New number of items to store.
         */
        virtual void resize(size_t new_size) = 0;

        /**
         * @brief Removes the item at @p index by swapping it with the last item and popping it (O(1)).
         * @param index Index of the item to remove.
         */
        virtual void erase_swap(size_t index) = 0; // Utile pour la suppression
    };

    /**
     * @class Variable
     * @brief Typed container attaching one value of type @p T to each cell of a given kind
     * (node, edge, face or 3D cell).
     * @tparam T Type of the values stored per cell.
     */
    /*----------------------------------------------------------------------------*/
    template<typename T>
    class Variable : public IVariable {
    public:
        /**
         * @brief Constructor
         * @param ASize initial size of the variable
         */
        explicit Variable(size_t ASize = 0) : m_data(ASize) {}

        /**
         * @brief Accesses the value attached to a cell by index (no bounds check).
         * @param AIdx Index of the cell.
         * @return Reference to the value attached to the cell.
         */
        T &operator[](size_t AIdx) { return m_data[AIdx]; }
        /**
         * @brief Accesses the value attached to a cell by index (no bounds check).
         * @param AIdx Index of the cell.
         * @return Const reference to the value attached to the cell.
         */
        const T &operator[](size_t AIdx) const { return m_data[AIdx]; }

        /**
         * @brief Accesses the value attached to a cell by index (bounds-checked).
         * @param idx Index of the cell.
         * @return Reference to the value attached to the cell.
         * @throw std::out_of_range if @p idx is out of bounds.
         */
        T &at(size_t idx) { return m_data.at(idx); }
        /**
         * @brief Accesses the value attached to a cell by index (bounds-checked).
         * @param idx Index of the cell.
         * @return Const reference to the value attached to the cell.
         * @throw std::out_of_range if @p idx is out of bounds.
         */
        const T &at(size_t idx) const { return m_data.at(idx); }

        // --- Functions required for range-based loops (for (auto& val : var)) ---
        /** @brief Iterator to the first stored value. @return Begin iterator. */
        auto begin() noexcept { return m_data.begin(); }
        /** @brief Iterator past the last stored value. @return End iterator. */
        auto end() noexcept { return m_data.end(); }
        /** @brief Const iterator to the first stored value. @return Begin const iterator. */
        auto begin() const noexcept { return m_data.begin(); }
        /** @brief Const iterator past the last stored value. @return End const iterator. */
        auto end() const noexcept { return m_data.end(); }

        /** @copydoc IVariable::resize */
        void resize(size_t new_size) override { m_data.resize(new_size); }
        /** @copydoc IVariable::size */
        size_t size() const override { return m_data.size(); }

        /** @copydoc IVariable::erase_swap */
        void erase_swap(size_t index) override {
            if (index < m_data.size()) {
                m_data[index] = std::move(m_data.back());
                m_data.pop_back();
            }
        }

    private:
        std::vector<T> m_data;
    };

    /**
     * @class VariableRegistry
     * @brief Owns and manages the named Variable<T> attached to a given kind of mesh cell
     * (node, edge, face or 3D cell), keeping all of them the same size.
     */
    /*----------------------------------------------------------------------------*/
    class VariableRegistry {
        // All the instanciation of  UnstructuredMesh<T> are friend
        template<typename TopologyTraits>
        friend class UnstructuredMesh;

    private:
        /**
         * @brief Constructor
         * @param ASize Initial size of the registry
         */
        explicit VariableRegistry(size_t ASize = 0) : m_current_size(ASize) {}

        /**
         * @brief Access to the variable of name @p AName, that contains items of type @p T.
         * @tparam T Type of the items to store.
         * @param AName Name given to the variable.
         * @return A reference onto the variable.
         * @throw An exception if such a variable does not exist.
         * @throw An exception if the variable named @p AName exists but is not of type @p T.
         */
        template<typename T>
        Variable<T> &get(const std::string &AName) {
            auto it = m_variables.find(AName);
            if (it == m_variables.end()) {
                throw std::out_of_range("Unfound variable: " + AName);
            }

            auto result = std::dynamic_pointer_cast<Variable<T>>(it->second);
            if (!result) {
                throw std::bad_cast(); //Wrong type required
            }
            return *result; // Returns a reference to a Variable<T> object
        }

        /**
         * @brief Access to the variable of name @p AName, that contains items of type @p T.
         * @tparam T Type of the items to store.
         * @param AName Name given to the variable.
         * @return A reference onto the variable.
         * @throw An exception if such a variable does not exist.
         * @throw An exception if the variable named @p AName exists but is not of type @p T.
        */
        template<typename T>
        const Variable<T> &get(const std::string &AName) const {
            auto it = m_variables.find(AName);
            if (it == m_variables.end()) {
                throw std::out_of_range("Attribut non trouvé : " + AName);
            }

            auto result = std::dynamic_pointer_cast<const Variable<T>>(it->second);
            if (!result) {
                throw std::bad_cast();
            }
            return *result;
        }

        /**
         * @brief Add a new variable of name @p AName containing items of type @p T.
         * @tparam T type of the stored items.
         * @param AName Name of the variable.
         * @return a reference onto the created variable.
         */
        template<typename T>
        Variable<T> &add(const std::string &AName) {
            auto ptr = std::make_shared<Variable<T>>(m_current_size);

            // Atomic insertion: check the existence and insert in one research
            auto [it, inserted] = m_variables.try_emplace(AName, ptr);

            if (!inserted) {
                throw std::invalid_argument("A variable of '" + AName + "' already exists in the registry.");
            }

            return *ptr;
        }

        /**
         * @brief Remove the variable of name @p AName.
         * @param AName Name of the variable.
         * @return true if the variable exists and is removed, false otherwise.
         */
        bool remove(const std::string &AName) { return m_variables.erase(AName) > 0; }

        /**
        * @brief Check if a variable of name @p AName is in the registry.
        * @param AName Name of the variable.
        * @return true if the variable exists, false otherwise.
        */
        bool has(const std::string &AName) const noexcept { return m_variables.contains(AName); }

        /**
         * @brief Clear the whole registry
         */
        void clear() { m_variables.clear(); }

        /**
         * @brief Resizes every registered variable to @p new_size (e.g. when growing the mesh).
         * @param new_size New number of items each variable must hold.
         */
        void resize_all(size_t new_size) {
            m_current_size = new_size;
            for (auto &[_, attr] : m_variables) {
                attr->resize(new_size); // Appelle IVariable::resize
            }
        }

        /**
         * @brief Removes the item at @p index from every registered variable (O(1) cell deletion),
         * by swapping it with the last item and popping it.
         * @param index Index of the item to remove.
         */
        void erase_swap_all(size_t index) {
            for (auto &[_, attr] : m_variables) {
                attr->erase_swap(index); // Appelle IVariable::erase_swap
            }
            --m_current_size;
        }

        /**
         * @brief Gets the number of items stored by each variable of this registry.
         * @return The current number of items.
         */
        size_t current_size() const noexcept { return m_current_size; }

        /**
         * @brief Gets the number of variables registered in this registry.
         * @return The number of registered variables.
         */
        size_t num_attributes() const noexcept { return m_variables.size(); }

    private:
        std::unordered_map<std::string, std::shared_ptr<IVariable>> m_variables;
        size_t m_current_size = 0;
    };

    /*----------------------------------------------------------------------------*/
} // namespace gecko
